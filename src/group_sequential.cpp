/*
 * Copyright © 2015-2019 Barcelona Supercomputing Center (BSC)
 *
 * This file is part of SMUFIN Core. SMUFIN Core is released under the SMUFIN
 * Public License, and may not be used except in compliance with it. This file
 * is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; see the SMUFIN Public License for more details. You should have
 * received a copy of the SMUFIN Public License along with this file. If not,
 * see <https://github.com/smufin/smufin-core/blob/master/COPYING>.
 *
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2018
 */

#include "group_sequential.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <rocksdb/db.h>

#include "db.hpp"
#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

group_sequential::group_sequential(const sm_config &conf) : stage(conf)
{
    init_mapping(conf, _conf.num_partitions, _conf.num_groupers,
                 _group_map_l1, _group_map_l2);
    _leads_size = _conf.leads_size / _conf.num_partitions / _conf.num_groupers;
    _executable["run"] = std::bind(&group_sequential::run, this);
    _executable["stats"] = std::bind(&group_sequential::stats, this);
}

void group_sequential::run()
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

    for (int i = 0; i < 2; i++) {
        _seq[i] = new seq_table();
        _k2i[i] = new k2i_table();
    }

    _seq[NN]->resize(420000000);
    _seq[TN]->resize(500000000);
    _k2i[NN]->resize(30000000);
    _k2i[TN]->resize(60000000);

    string dir;
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
        p = decode_pos(it->value().ToString());

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
        map_mer(m);
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

    // 2. Iterate K2I retrieving all k-mers seen in candidate positions.

    std::vector<sm_idx_set> sets = {NN, TN};
    for (auto& set: sets) {
        start = std::chrono::system_clock::now();
        rdb_handle k2i;
        open_index_full_iter(_conf, K2I, set, k2i);

        int num_seen = 0;
        int num_kmer = 0;
        istart = std::chrono::system_clock::now();
        it = k2i.db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) {
            num_kmer++;

            string kmer = it->key().ToString();
            string list = it->value().ToString();
            k2i_table::const_iterator kit = _k2i[set]->find(kmer);
            if (kit != _k2i[set]->end()) {
                (*_k2i[set])[kmer] = list;
                num_seen++;
            }

            if (num_kmer % 1000000 == 0) {
                iend = std::chrono::system_clock::now();
                itime = iend - istart;
                cout << "K: " << sm::sets[set] << " " << num_kmer << " "
                     << num_seen << " " << itime.count() << endl;
                istart = std::chrono::system_clock::now();
            }
        }

        cout << "Number of kmers (" << sm::sets[set] << "): "
             << _k2i[set]->size() << endl;
        cout << "Number of kmers seen (" << sm::sets[set] << "): "
             << num_seen << endl;
        cout << "Number of IDs seen (" << sm::sets[set] << "): "
             << _seq[set]->size() << endl;


        delete it;
        delete k2i.cfs[0];
        delete k2i.db;

        end = std::chrono::system_clock::now();
        time = end - start;
        cout << "Kmer iteration time (" << sm::sets[set] << "): "
             << time.count() << endl;
    }

    // 3. Iterate & collect reads.

    for (auto& set: sets) {
        start = std::chrono::system_clock::now();

        rdb_handle seq;
        open_index_full_iter(_conf, SEQ, set, seq);

        int num_seen = 0;
        int num_read = 0;
        istart = std::chrono::system_clock::now();
        it = seq.db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) {
            num_read++;

            string sid = it->key().ToString();
            string read_str = it->value().ToString();
            sm_read_code read;
            encode_read(read_str, read);
            (*_seq[set])[sid] = read;

            if (num_read % 10000000 == 0) {
                iend = std::chrono::system_clock::now();
                itime = iend - istart;
                cout << "R: " << sm::sets[set] << " " << num_read << " "
                     << num_seen << " " << itime.count() << endl;
                istart = std::chrono::system_clock::now();
            }
        }

        delete it;
        delete seq.cfs[0];
        delete seq.db;

        end = std::chrono::system_clock::now();
        time = end - start;
        cout << "Read iteration time (" << sm::sets[set] << "): "
             << time.count() << endl;
    }

    // 4. Populate candidate groups.

    spawn("populate", std::bind(&group_sequential::populate, this,
          std::placeholders::_1), _conf.num_groupers);
}

void group_sequential::encode_read(std::string& str, sm_read_code& read)
{
    read.len = str.size();
    for (int i = 0; i < read.len; i += 32) {
        read.seq[i/32] = strtob4(str.substr(i, 32).c_str());
    }
}

void group_sequential::decode_read(sm_read_code& read, std::string& str)
{
    char decode[4] = {'A', 'C', 'G', 'T'};
    str = "";

    int len_bits = read.len * 2;
    int pos = len_bits / 64;
    int shift = 64 - (len_bits % 64);
    read.seq[pos] = read.seq[pos] << shift;

    int l = 0;
    int b = 0;
    while (l < read.len) {
        int pos = b / 64;
        int shift = 62 - (b % 64);
        int i = ((read.seq[pos] >> shift) & 3);
        char c = decode[i];
        str += c;
        l++;
        b += 2;
    }
}

void group_sequential::select_candidate(int gid, string& sid, string& seq,
                                        string& dseq, std::vector<int>& pos,
                                        int dir)
{
    std::vector<string> kmers;
    for (int p: pos) {
        string kmer = dseq.substr(p, _conf.k);
        kmers.push_back(kmer);
        k2i_table::const_iterator it = _k2i[0]->find(kmer);
        if (it == _k2i[0]->end()) {
            (*_k2i[0])[kmer] = string();
            (*_k2i[1])[kmer] = string();
        }
    }
    (*_l2r[gid])[sid] = seq;
    (*_l2p[gid])[sid][dir] = pos;
    (*_l2k[gid])[sid][dir] = kmers;
}

void group_sequential::populate(int gid)
{
    const char comp_code[] = "ab";
    const char kind_code[] = "nt";

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    std::ofstream ofs;
    string file = _conf.output_path_group + "/group." +
                  std::to_string(_conf.pid) + "-" + std::to_string(gid) +
                  ".json";
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

        populate_index(gid, lid, kmers[0], NN, keep, drop);
        populate_index(gid, lid, kmers[0], TN, keep, drop);
        populate_index(gid, lid, kmers[1], NN, keep, drop);
        populate_index(gid, lid, kmers[1], TN, keep, drop);

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
                seq_table::const_iterator sit = _seq[i]->find(sid);
                if (sit == _seq[i]->end())
                    continue;
                if (!first_read)
                    ofs << ",";
                first_read = false;
                sm_read_code read = sit->second;
                seq = "";
                decode_read(read, seq);
                ofs << "[\"" << sid << "\",\"" << seq << "\"]";
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

void group_sequential::populate_index(int gid, const string& lid,
                                      const std::vector<string>& kmers,
                                      int kind, kmer_count& keep,
                                      kmer_count& drop)
{
    for (string kmer: kmers) {
        k2i_table::const_iterator it = _k2i[kind]->find(kmer);
        if (it == _k2i[kind]->end()) {
            continue;
        }

        string list = it->second;
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

void group_sequential::stats()
{
    uint64_t total = 0;

    for (int i = 0; i < _conf.num_groupers; i++) {
        cout << "Groups " << i << ": " << _num_groups[i] << endl;
        total += _num_groups[i];
    }

    cout << "Number of groups: " << total << endl;
}
