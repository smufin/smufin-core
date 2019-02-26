/*
 * Copyright © 2015-2018 Barcelona Supercomputing Center (BSC)
 *
 * This file is part of SMUFIN Core. SMUFIN Core is released under the SMUFIN
 * Public License, and may not be used except in compliance with it. This file
 * is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; see the SMUFIN Public License for more details. You should have
 * received a copy of the SMUFIN Public License along with this file. If not,
 * see <https://github.com/smufin/smufin-core/blob/master/COPYING>.
 *
 * Jordà Polo <jorda.polo@bsc.es>, 2018
 */

#include "group_rocks.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "db.hpp"
#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

group_rocks::group_rocks(const sm_config &conf) : stage(conf)
{
    init_mapping(conf, _conf.num_partitions, _conf.num_groupers,
                 _group_map_l1, _group_map_l2);
    _leads_size = _conf.leads_size / _conf.num_partitions / _conf.num_groupers;
    _executable["run"] = std::bind(&group_rocks::run, this);
    _executable["stats"] = std::bind(&group_rocks::stats, this);
}

void group_rocks::run()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    for (int i = 0; i < _conf.num_groupers; i++) {
        open_groups_part(_conf, i, _groups[i]);
    }

    string dir;
    rocksdb::Status status;

    rdb_handle i2p_tm;
    rdb_handle seq_tm;
    open_index_full_iter(_conf, I2P, TM, i2p_tm);
    open_index_full_iter(_conf, SEQ, TM, seq_tm);

    const int min = _conf.window_min;
    const int len = _conf.window_len;

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

        sm_group group;
        bool match = false;

        if (a_len >= KMIN && a_len <= KMAX && match_window(a_pos, min, len)) {
            select_candidate(gid, sid, read, read, a_pos, 0, group);
            num_match_a++;
            match = true;
        }

        if (b_len >= KMIN && b_len <= KMAX && match_window(b_pos, min, len)) {
            char buf[MAX_READ_LEN + 1];
            copy(read.begin(), read.end(), buf);
            buf[read_length] = '\0';
            revcomp(buf, read_length);
            string directed_read = string(buf);
            select_candidate(gid, sid, read, directed_read, b_pos, 1, group);
            num_match_b++;
            match = true;
        }

        if (match) {
            std::stringstream buf;
            msgpack::pack(buf, group);
            rocksdb::WriteOptions w_options;
            w_options.disableWAL = true;
            _groups[gid].db->Put(w_options, _groups[gid].cfs[0], sid, buf.str());
            _num_groups[gid]++;
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
        uint64_t n = 0;
        _groups[i].db->GetAggregatedIntProperty("rocksdb.estimate-num-keys", &n);
        cout << "Number of candidates (" << std::to_string(i) << "): "
             << _num_groups[i] << " " << n << endl;
    }

    // 2. Populate candidate groups with kmers and reads.

    open_index_full_read(_conf, K2I, NN, _k2i[NN]);
    open_index_full_read(_conf, K2I, TN, _k2i[TN]);
    open_index_full_read(_conf, SEQ, NN, _seq[NN]);
    open_index_full_read(_conf, SEQ, TN, _seq[TN]);

    spawn("populate", std::bind(&group_rocks::populate, this,
          std::placeholders::_1), _conf.num_groupers);

    delete _k2i[NN].cfs[0], _k2i[TN].cfs[0], _seq[NN].cfs[0], _seq[TN].cfs[0];
    delete _k2i[NN].db, _k2i[TN].db, _seq[NN].db, _seq[TN].db;

    // 3. Write results.

    spawn("dump", std::bind(&group_rocks::dump, this, std::placeholders::_1),
          _conf.num_groupers);
}

void group_rocks::select_candidate(int gid, string& sid, string& seq,
                                   string& dseq, std::vector<int>& pos,
                                   int dir, sm_group& group)
{
    group.lead = sm_group_read(sid, seq);

    std::vector<int> pos_v;
    std::vector<sm_group_kmer> kmers_v;

    for (int p: pos) {
        pos_v.push_back(p);
        string kmer = dseq.substr(p, _conf.k);
        kmers_v.push_back(sm_group_kmer(kmer, std::array<int, 4>()));
    }

    group.pos[dir] = pos_v;
    group.kmers[dir] = kmers_v;
}

void group_rocks::populate(int gid)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    rocksdb::Status status;
    rocksdb::ReadOptions r_options;
    rocksdb::WriteOptions w_options;
    w_options.disableWAL = true;

    uint64_t num_groups = 0;

    rocksdb::Iterator* it;
    it = _groups[gid].db->NewIterator(r_options, _groups[gid].cfs[0]);
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        num_groups++;

        string val = it->value().ToString();
        msgpack::object_handle oh = msgpack::unpack(val.c_str(), val.size());
        msgpack::object obj = oh.get();
        sm_group group;
        obj.convert(group);

        populate_kmers(group, NN, _k2i[NN]);
        populate_kmers(group, TN, _k2i[TN]);

        populate_reads(group, NN, _seq[NN]);
        populate_reads(group, TN, _seq[TN]);

        std::stringstream buf;
        msgpack::pack(buf, group);
        _groups[gid].db->Put(w_options, _groups[gid].cfs[0], group.lead.first,
                             buf.str());

        if (num_groups % 100 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "P: " << gid << " " << num_groups << " " << time.count()
                 << "\n";
            start = std::chrono::system_clock::now();
        }
    }

    delete it;
}

void group_rocks::populate_kmers(sm_group& group, sm_idx_set set,
                                 rdb_handle &rdb)
{
    rocksdb::Status status;
    std::vector<sm_dir> dirs= {DIR_A, DIR_B};
    for (auto& dir: dirs) {
        for (auto& k: group.kmers[dir]) {
            string kmer = k.first;
            string list;
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
                k.second[2 + set] += sids.size();
                continue;
            }

            k.second[set] += sids.size();
            for (const auto& sid: sids) {
                group.reads[set].push_back(std::pair<string,string>(sid, ""));
            }
        }
    }
}

void group_rocks::populate_reads(sm_group& group, sm_idx_set set,
                                 rdb_handle &rdb)
{
    rocksdb::Status status;
    std::set<sm_group_read> sids(group.reads[set].begin(),
                                 group.reads[set].end());
    std::vector<sm_group_read> reads;
    for (auto& k: sids) {
        string sid = k.first;
        string seq;
        status = rdb.db->Get(rocksdb::ReadOptions(), sid, &seq);
        if (!status.ok()) {
            continue;
        }
        reads.push_back(std::pair<string, string>(sid, seq));
    }
    group.reads[set].clear();
    group.reads[set] = reads;
}

void group_rocks::dump(int gid)
{
    std::ofstream ofs;
    string file = _conf.output_path_group + "/group." +
                  std::to_string(_conf.pid) + "-" + std::to_string(gid) +
                  ".msgpack";
    ofs.open(file);

    rocksdb::ReadOptions r_options;
    rocksdb::Iterator* it;
    it = _groups[gid].db->NewIterator(r_options, _groups[gid].cfs[0]);
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        string val = it->value().ToString();
        ofs << val;
    }
    delete it;

    ofs.close();
}

void group_rocks::stats()
{
    uint64_t total = 0;

    for (int i = 0; i < _conf.num_groupers; i++) {
        cout << "Groups " << i << ": " << _num_groups[i] << endl;
        total += _num_groups[i];
    }

    cout << "Number of groups: " << total << endl;
}
