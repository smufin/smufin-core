#include "group.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

group::group(const sm_config &conf) : stage(conf)
{
    init_mapping(conf, _conf.num_partitions, _conf.num_groupers,
                 _group_map_l1, _group_map_l2);
    _leads_size = _conf.leads_size / _conf.num_partitions / _conf.num_groupers;
    _executable["run"] = std::bind(&group::run, this);
    _executable["stats"] = std::bind(&group::stats, this);
}

void group::run()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    for (int i = 0; i < _conf.num_groupers; i++) {
        _l2p[i] = new l2p_table(_leads_size);
        _l2k[i] = new l2k_table(_leads_size);
        _l2i[i] = new l2i_table(_leads_size);
        _l2r[i] = new l2r_table(_leads_size);
    }

    rocksdb::Status status;

    rdb_handle i2p_tm;
    rdb_handle seq_tm;
    open_index_full_iter(_conf, I2P, TM, i2p_tm);
    open_index_full_iter(_conf, SEQ, TM, seq_tm);

    int min = _conf.window_min;
    int len = _conf.window_len;

    string sid;
    sm_pos_bitmap p;

    int num_all = 0;
    int num_exist = 0;
    int num_map = 0;
    int num_match_a = 0;
    int num_match_b = 0;

    std::chrono::time_point<std::chrono::system_clock> istart, iend;
    std::chrono::duration<double> itime;
    istart = std::chrono::system_clock::now();

    // 1. Iterate through candidates.

    rocksdb::Iterator* it = i2p_tm.db->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        num_all++;

        sid = it->key().ToString();
        p = decode_pos(it->value().data());

        std::vector<int> a_pos;
        std::vector<int> b_pos;

        string read;
        rocksdb::Status status;
        status = seq_tm.db->Get(rocksdb::ReadOptions(), sid, &read);
        if (!status.ok())
            continue;

        if (num_all % 100000 == 0) {
            iend = std::chrono::system_clock::now();
            itime = iend - istart;
            cout << "S: " << num_all << " " << itime.count() << endl;
            istart = std::chrono::system_clock::now();
        }

        num_exist++;

        int read_length = read.size();

        if (read_length == 0)
            continue;

        if (read.find("N") != std::string::npos)
            continue;

        uint64_t m = 0;
        string sub = read.substr(0, MAP_LEN);
        memcpy(&m, sub.c_str(), MAP_LEN);
        hash_5mer(m);
        if (_group_map_l1[m] != _conf.pid)
            continue;
        int gid = _group_map_l2[m];

        num_map++;

        get_positions(p.a, &a_pos);
        get_positions(p.b, &b_pos);

        int a_len = a_pos.size();
        int b_len = b_pos.size();

        if (a_len >= KMIN && a_len <= KMAX && match_window(a_pos, min, len)) {
            select_candidate(gid, sid, read, read, a_pos, 0);
            num_match_a++;
        }

        if (b_len >= KMIN && b_len <= KMAX && match_window(b_pos, min, len)) {
            char buf[MAX_READ_LEN + 1];
            copy(read.begin(), read.end(), buf);
            buf[read_length] = '\0';
            revcomp(buf, read_length);
            string directed_read = string(buf);
            select_candidate(gid, sid, read, directed_read, b_pos, 1);
            num_match_b++;
        }
    }

    delete it;
    delete i2p_tm.cfs[0], seq_tm.cfs[0];
    delete i2p_tm.db, seq_tm.db;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Candidate selection time: " << time.count() << endl;

    cout << "Number of iterated I2P: " << num_all << endl;
    cout << "Number of existing I2P: " << num_exist << endl;
    cout << "Number of mapped I2P: " << num_map << endl;
    cout << "Number of A-matching I2P: " << num_match_a << endl;
    cout << "Number of B-matching I2P: " << num_match_b << endl;
    for (int i = 0; i < _conf.num_groupers; i++) {
        cout << "Number of candidates (" << std::to_string(i) << "): "
             << _l2k[i]->size() << endl;
    }

    // 2. Populate candidate groups.

    spawn("populate", std::bind(&group::populate, this, std::placeholders::_1),
          _conf.num_groupers);
}

void group::select_candidate(int gid, string& sid, string& seq, string& dseq,
                             std::vector<int>& pos, int dir)
{
    std::vector<string> kmers;
    for (int p: pos) {
        if (p > 50)
            break;
        string kmer = dseq.substr(p, _conf.k);
        kmers.push_back(kmer);
    }
    (*_l2r[gid])[sid] = seq;
    (*_l2p[gid])[sid][dir] = pos;
    (*_l2k[gid])[sid][dir] = kmers;
}

void group::populate(int gid)
{
    const char comp_code[] = "ab";
    const char kind_code[] = "nt";

    rdb_handle _seq[2];
    rdb_handle _k2i[2];
    open_index_full_read(_conf, K2I, NN, _k2i[NN]);
    open_index_full_read(_conf, K2I, TN, _k2i[TN]);
    open_index_full_read(_conf, SEQ, NN, _seq[NN]);
    open_index_full_read(_conf, SEQ, TN, _seq[TN]);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    rocksdb::Status status;

    std::ofstream ofs;
    string file = _conf.output_path + "/group." + std::to_string(_conf.pid) +
                  "-" + std::to_string(gid) + ".json";
    ofs.open(file);
    ofs << "{";

    uint64_t num_groups = 0;
    bool first_group = true;
    for (const auto& it: *_l2k[gid]) {
        const string& lid = it.first;
        const k_value& kmers = it.second;

        kmer_count keep;
        kmer_count drop;

        for (int i = 0; i < 2; i++) {
            for (string kmer: kmers[i]) {
                keep[i][kmer] = 0;
                drop[i][kmer] = 0;
            }
        }

        populate_index(gid, lid, kmers[0], NN, keep, drop, _k2i[NN]);
        populate_index(gid, lid, kmers[0], TN, keep, drop, _k2i[TN]);
        populate_index(gid, lid, kmers[1], NN, keep, drop, _k2i[NN]);
        populate_index(gid, lid, kmers[1], TN, keep, drop, _k2i[TN]);

        l2r_table::const_iterator lit = _l2r[gid]->find(lid);
        if (lit == _l2r[gid]->end())
            continue;
        string seq = lit->second;

        if (!first_group)
            ofs << ",";
        first_group = false;

        ofs << "\"" << lid << "\":{";
        ofs  << "\"lead\":["
             << "\"" << lid << "\","
             << "\"" << seq << "\""
             << "],";

        for (int i = 0; i < 2; i++) {
            ofs << "\"pos-" << comp_code[i] << "\":[";
            bool first_pos = true;
            for (int p: (*_l2p[gid])[lid][i]) {
                if (!first_pos)
                    ofs << ",";
                first_pos = false;
                ofs << p;
            }
            ofs << "],";

            ofs << "\"kmers-" << comp_code[i] << "\":[";
            bool first_kmer = true;
            for (string kmer: kmers[i]) {
                int kept_n = keep[0][kmer];
                int dropped_n = drop[0][kmer];
                int kept_t = keep[1][kmer];
                int dropped_t = drop[1][kmer];
                if (!first_kmer)
                    ofs << ",";
                first_kmer = false;
                ofs << "[\"" << kmer << "\"," << kept_n << "," << kept_t << ","
                     << dropped_n << "," << dropped_t << "]";
            }
            ofs << "],";
        }

        for (int i = 0; i < 2; i++) {
            bool first_read = true;
            ofs << "\"reads-" << kind_code[i] << "\":[";
            for (string sid: (*_l2i[gid])[lid][i]) {
                string read;
                status = _seq[i].db->Get(rocksdb::ReadOptions(), sid, &read);
                if (!status.ok())
                    continue;
                if (!first_read)
                    ofs << ",";
                first_read = false;
                ofs << "[\"" << sid << "\",\"" << read << "\"]";
            }
            ofs << ( i == 1 ? "]" : "]," );
        }

        ofs << "}";
        num_groups++;

        if (num_groups % 100 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "P: " << gid << " " << num_groups << " "
                 << time.count() << "\n";
            start = std::chrono::system_clock::now();
        }
    }

    ofs << "}";
    _num_groups[gid] = num_groups;
}

void group::populate_index(int gid, const string& lid,
                           const std::vector<string>& kmers, int kind,
                           kmer_count& keep, kmer_count& drop,
                           rdb_handle &rdb)
{
    for (string kmer: kmers) {
        string list;
        rocksdb::Status status;
        status = rdb.db->Get(rocksdb::ReadOptions(), kmer, &list);
        if (!status.ok()) {
            continue;
        }

        boost::trim_if(list, boost::is_any_of(" "));
        if (list.size() == 0) {
            continue;
        }

        std::unordered_set<string> sids;
        boost::split(sids, list, boost::is_any_of(" "));

        if (sids.size() > _conf.max_group_reads) {
            drop[kind][kmer] += sids.size();
            continue;
        }

        (*_l2i[gid])[lid][kind].insert(sids.begin(), sids.end());
        keep[kind][kmer] += sids.size();
    }
}

void group::stats()
{
    uint64_t total = 0;

    for (int i = 0; i < _conf.num_groupers; i++) {
        cout << "Groups " << i << ": " << _num_groups[i] << endl;
        total += _num_groups[i];
    }

    cout << "Number of groups: " << total << endl;
}

void get_positions(const uint64_t bitmap[POS_LEN], std::vector<int> *pos)
{
    for (int i = 0; i < 2; i++) {
        unsigned long tmp = bitmap[i];
        int offset = i * 64;
        while (tmp > 0) {
            int p = __builtin_ffsl(tmp) - 1;
            tmp &= (tmp - 1);
            pos->push_back(p + offset);
        }
    }
}

bool match_window(const std::vector<int> pos, int window_min, int window_len)
{
    if (pos.size() < window_min) {
        return false;
    }

    for (int i = 0; i <= pos.size() - window_min; i++) {
        std::vector<int> sub(pos.begin() + i, pos.begin() + i + window_min);
        if (sub.back() - sub.front() < window_len) {
            return true;
        }
    }

    return false;
}
