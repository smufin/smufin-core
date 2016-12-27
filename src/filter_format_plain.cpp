#include "filter_format_plain.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

using std::cout;
using std::endl;
using std::string;

void filter_format_plain::update(kseq_t *seq, int pos, bool rev, char kmer[],
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

bool filter_format_plain::flush()
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

void filter_format_plain::stats()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    cout << "Size SEQ: " << _ids[NN].size() << " " << _ids[TN].size() << " "
         << _ids[TM].size() << endl;
    cout << "Size K2I: " << _k2i[NN].size() << " " << _k2i[TN].size() << " "
         << _k2i[TM].size() << endl;
    cout << "Size I2P: " << _i2p[TM].size() << endl;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Time filter/stats: " << time.count() << endl;
}

void filter_format_plain::dump()
{
    std::vector<std::thread> seqs;
    for (int i = 0; i < NUM_SETS; i++)
        seqs.push_back(std::thread(&filter_format_plain::write_seq, this, i));
    cout << "Spawned " << seqs.size() << " SEQ writer threads" << endl;

    std::vector<std::thread> k2is;
    for (int i = 0; i < NUM_SETS; i++)
        k2is.push_back(std::thread(&filter_format_plain::write_k2i, this, i));
    cout << "Spawned " << k2is.size() << " K2I writer threads" << endl;

    std::thread i2p = std::thread(&filter_format_plain::write_i2p, this, TM);
    cout << "Spawned I2P writer thread" << endl;
    i2p.join();

    for (auto& k2i: k2is)
        k2i.join();
    for (auto& seq: seqs)
        seq.join();
}

void filter_format_plain::write_seq(int set)
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

void filter_format_plain::write_k2i(int set)
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

void filter_format_plain::write_i2p(int set)
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
