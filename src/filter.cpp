#include "filter.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <thread>

#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::string;

filter::filter(const sm_config &conf) : stage(conf)
{
    std::ifstream ifs(_conf.input_file);
    if (!ifs.good()) {
        cout << "Invalid input file" << endl;
        exit(1);
    }

    for (string line; std::getline(ifs, line);) {
        _input_queue.enqueue(line);
        _input_count++;
    }

    _format = new filter_format(conf);

    _executable["run"] = std::bind(&filter::run, this);
    _executable["dump"] = std::bind(&filter::dump, this);
    _executable["stats"] = std::bind(&filter::stats, this);
}

void filter::chain(const stage* prev)
{
    cout << "Chain: filter" << endl;
    _count = static_cast<const count*>(prev);
}

void filter::run()
{
    spawn("filter", std::bind(&filter::load, this, std::placeholders::_1),
          _conf.num_filters);
}

void filter::stats()
{
    _format->stats();
}

void filter::dump()
{
    _format->dump();
}

void filter::load(int fid)
{
    string file;
    while (_input_count > 0) {
        while (_input_queue.try_dequeue(file)) {
            load_file(fid, file);
            _input_count--;
        }
    }
}

void filter::load_file(int fid, string file)
{
    // Identify read kind from file name.
    sm_read_kind kind = NORMAL_READ;
    std::size_t found = file.find("_T_");
    if (found != std::string::npos)
        kind = CANCER_READ;

    int len;
    int nreads = 0;
    gzFile in = gzopen(file.c_str(), "rb");

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    kseq_t *seq = kseq_init(in);
    while ((len = kseq_read(seq)) >= 0) {
        nreads++;

        if (lq_count(seq->qual.s, seq->qual.l) > len/10)
            continue;

        int p = 0;
        int l = seq->seq.l;
        int n = seq->seq.l;
        char *ps;

        while ((ps = (char*) memchr(&seq->seq.s[p], 'N', l - p)) != NULL) {
            n = ps - &seq->seq.s[p];
            if (n > 0) {
                if (kind == CANCER_READ)
                    filter_cancer(fid, seq, &seq->seq.s[p], n);
                else
                    filter_normal(fid, seq, &seq->seq.s[p], n);
                p += n;
            }
            p++;
        }

        n = l - p;
        if (n > 0) {
            if (kind == CANCER_READ)
                filter_cancer(fid, seq, &seq->seq.s[p], n);
            else
                filter_normal(fid, seq, &seq->seq.s[p], n);
        }

        if (nreads % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "F: " << fid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }

        if (nreads % 10000000 == 0) {
            bool f = _format->flush();
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "W: " << fid << " " << time.count() << " " << f << endl;
            start = std::chrono::system_clock::now();
        }
    }

    kseq_destroy(seq);
    gzclose(in);
}

void filter::filter_normal(int fid, kseq_t *seq, const char *sub, int len)
{
    if (len < KMER_LEN)
        return;

    char kmer[KMER_LEN + 1];
    for (int i = 0; i <= len - KMER_LEN; i++) {
        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';
        filter_all(fid, seq, i, false, kmer, NN);

        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';
        krevcomp(kmer);
        filter_all(fid, seq, i, true, kmer, NN);
    }
}

void filter::filter_cancer(int fid, kseq_t *seq, const char *sub, int len)
{
    if (len < KMER_LEN)
        return;

    char kmer[KMER_LEN + 1];
    for (int i = 0; i <= len - KMER_LEN; i++) {
        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';
        filter_branch(fid, seq, i, false, kmer, TM);
        filter_all(fid, seq, i, false, kmer, TN);

        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';
        krevcomp(kmer);
        filter_branch(fid, seq, i, true, kmer, TM);
        filter_all(fid, seq, i, true, kmer, TN);
    }
}

int filter::get_value(int fid, char kmer[], sm_table::const_iterator *it)
{
    char last = kmer[KMER_LEN - 1];
    kmer[KMER_LEN - 1] = '\0';

    uint64_t m = 0;
    memcpy(&m, &kmer[1], MAP_LEN);
    hash_5c_map(m);

    if (map_l1[m] != _conf.pid)
        return -1;
    int sid = map_l2[m];
    sm_key key = strtob4(&kmer[1]);

    const sm_table* table = (*_count)[sid];
    *it = table->find(key);
    if (*it == table->end())
        return -1;

    kmer[KMER_LEN - 1] = last;
    return 0;
}

void filter::filter_branch(int fid, kseq_t *seq, int pos, bool rev,
                           char kmer[], sm_set set)
{
    sm_table::const_iterator it;
    if (get_value(fid, kmer, &it) != 0)
        return;
    int f = sm::code[kmer[0]] - '0';
    int l = sm::code[kmer[KMER_LEN - 1]] - '0';
    uint32_t nc = it->second.v[f][l][NORMAL_READ];
    uint32_t tc = it->second.v[f][l][CANCER_READ];
    uint32_t nsum = 0;
    uint32_t tsum = 0;
    for (l = 0; l < 4; l++) {
        nsum += it->second.v[f][l][NORMAL_READ];
        tsum += it->second.v[f][l][CANCER_READ];
    }
    filter_kmer(seq, pos, rev, kmer, nc, tc, nsum, tsum, set);
}

void filter::filter_all(int fid, kseq_t *seq, int pos, bool rev,
                        char kmer[], sm_set set)
{
    sm_table::const_iterator it;
    if (get_value(fid, kmer, &it) != 0)
        return;
    for (int f = 0; f < 4; f++) {
        kmer[0] = sm::alpha[f];
        uint32_t nsum = 0;
        uint32_t tsum = 0;
        for (int l = 0; l < 4; l++) {
            nsum += it->second.v[f][l][NORMAL_READ];
            tsum += it->second.v[f][l][CANCER_READ];
        }
        for (int l = 0; l < 4; l++) {
            kmer[KMER_LEN - 1] = sm::alpha[l];
            uint32_t nc = it->second.v[f][l][NORMAL_READ];
            uint32_t tc = it->second.v[f][l][CANCER_READ];
            filter_kmer(seq, pos, rev, kmer, nc, tc, nsum, tsum, set);
        }
    }
}

void filter::filter_kmer(kseq_t *seq, int pos, bool rev, char kmer[],
                         uint32_t nc, uint32_t tc, uint32_t nsum,
                         uint32_t tsum, sm_set set)
{
    if (tc >= _conf.min_tc && nc <= _conf.max_nc) {
        _format->update(seq, pos, rev, kmer, set);
    }
}

void filter_format::update(kseq_t *seq, int pos, bool rev, char kmer[],
                           sm_set set)
{
    char buf[512] = {0};
    sprintf(buf, "%s %s", seq->name.s, seq->seq.s);

    _mutex[set].lock();
    auto result = _ids[set].insert(seq->name.s);
    if (result.second == true) {
        _seq[set].insert(buf);
    }
    if (set == TM) {
        if (!rev)
            _i2p[set][seq->name.s].a[pos / 64] |= 1UL << (pos % 64);
        else
            _i2p[set][seq->name.s].b[pos / 64] |= 1UL << (pos % 64);
    }
    if (_k2i[set][kmer].size() <= _conf.max_k2i_reads) {
        _k2i[set][kmer].insert(seq->name.s);
    }
    _mutex[set].unlock();
}

bool filter_format::flush()
{
    bool flushed = false;
    for (int i = 0; i < NUM_SETS; i++) {
        _mutex[i].lock();
        if (_seq[i].size() > 1000000) {
            write_seq(i);
            _seq[i].clear();
            _seq[i] = std::unordered_set<std::string>();
            flushed = true;
        }
        _mutex[i].unlock();
    }
    return flushed;
}

void filter_format::stats()
{
    cout << "Filtered NN reads: " << _ids[NN].size() << endl;
    cout << "Filtered TN reads: " << _ids[TN].size() << endl;
    cout << "Filtered TM reads: " << _ids[TM].size() << endl;
}

void filter_format::dump()
{
    std::vector<std::thread> seqs;
    for (int i = 0; i < NUM_SETS; i++)
        seqs.push_back(std::thread(&filter_format::write_seq, this, i));
    cout << "Spawned " << seqs.size() << " SEQ writer threads" << endl;

    std::vector<std::thread> k2is;
    for (int i = 0; i < NUM_SETS; i++)
        k2is.push_back(std::thread(&filter_format::write_k2i, this, i));
    cout << "Spawned " << k2is.size() << " K2I writer threads" << endl;

    std::thread i2p = std::thread(&filter_format::write_i2p, this, TM);
    cout << "Spawned I2P writer thread" << endl;
    i2p.join();

    for (auto& k2i: k2is)
        k2i.join();
    for (auto& seq: seqs)
        seq.join();
}

void filter_format::write_seq(int set)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path << "/filter-seq-" << sm::sets[set] << "."
         << _conf.pid << ".txt";
    ofs.open(file.str(), std::ofstream::app);
    for (auto const &s: _seq[set]) {
        ofs << s << "\n";
    }
    ofs.close();
}

void filter_format::write_k2i(int set)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path << "/filter-k2i-" << sm::sets[set] << "."
         << _conf.pid << ".txt";
    ofs.open(file.str());
    for (auto const &kv: _k2i[set]) {
        if (kv.second.size() > _conf.max_k2i_reads)
            continue;
        ofs << kv.first << " " << kv.second.size();
        for (auto const &sid: kv.second) {
            ofs << " " << sid;
        }
        ofs << "\n";
    }
    ofs.close();
}

void filter_format::write_i2p(int set)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path << "/filter-i2p-" << sm::sets[set] << "."
         << _conf.pid << ".txt";
    ofs.open(file.str());
    for (auto const &kv: _i2p[set]) {
        const sm_pos_bitmap *p = &kv.second;
        ofs << kv.first << " ";
        ofs << p->a[0] << " " << p->a[1] << " ";
        ofs << p->b[0] << " " << p->b[1] << "\n";
    }
    ofs.close();
}
