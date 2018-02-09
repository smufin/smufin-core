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
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2018
 */

#include "index_format_rocks.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

#include "db.hpp"
#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

index_format_rocks::index_format_rocks(const sm_config &conf)
    : index_format(conf)
{
    int pid = _conf.pid;
    for (int iid = 0; iid < _conf.num_indexes; iid++) {
        for (auto set: {NN, TN, TM})
            open_index_part_load(_conf, SEQ, set, pid, iid, _seq[set][iid]);

        for (auto set: {NN, TN})
            open_index_part_load(_conf, K2I, set, pid, iid, _k2i[set][iid]);

        open_index_part_load(_conf, I2P, TM, pid, iid, _i2p[iid]);
    }
}

void index_format_rocks::update(int fid, const sm_read *read, int pos,
                                char kmer[], sm_dir dir, sm_idx_set set)
{
    string sid = read->id;
    rocksdb::WriteOptions options;
    options.disableWAL = true;
    int iid = fid % _conf.num_indexes;

    _seq[set][iid].db->Put(options, _seq[set][iid].cfs[0], sid, read->seq);

    if (set == TM) {
        sm_pos_bitmap p;
        string serialized;

        if (dir == DIR_A)
            p.a[pos / 64] |= 1UL << (pos % 64);
        else
            p.b[pos / 64] |= 1UL << (pos % 64);

        encode_pos(p, serialized);
        _i2p[iid].db->Merge(options, _i2p[iid].cfs[0], sid, serialized);
    } else {
        // TODO: Honour _conf.max_filter_reads
        std::stringstream ss;
        ss << sid << " ";
        _k2i[set][iid].db->Merge(options, _k2i[set][iid].cfs[0], kmer, ss.str());
    }
}

void index_format_rocks::dump()
{
    std::vector<rocksdb::DB*> list;
    for (int iid = 0; iid < _conf.num_indexes; iid++) {
        for (auto set: {NN, TN, TM})
            list.push_back(_seq[set][iid].db);
        for (auto set: {NN, TN})
            list.push_back(_k2i[set][iid].db);
        list.push_back(_i2p[iid].db);
    }
    spawn<rocksdb::DB*>("compact", std::bind(&index_format_rocks::compact,
                        this, std::placeholders::_1), list);
}

void index_format_rocks::compact(rocksdb::DB* db)
{
    db->CompactRange(rocksdb::CompactRangeOptions(), nullptr, nullptr);
}

void index_format_rocks::stats()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    rocksdb::Iterator* it;
    uint64_t seq_nn, seq_tn, seq_tm;
    uint64_t k2i_nn, k2i_tn;
    uint64_t i2p_tm;

    seq_nn = seq_tn = seq_tm = k2i_nn = k2i_tn = i2p_tm = 0;

    for (int iid = 0; iid < _conf.num_indexes; iid++) {
        it = _seq[NN][iid].db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) seq_nn++;
        it = _seq[TN][iid].db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) seq_tn++;
        it = _seq[TM][iid].db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) seq_tm++;

        it = _k2i[NN][iid].db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) k2i_nn++;
        it = _k2i[TN][iid].db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) k2i_tn++;

        it = _i2p[iid].db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) i2p_tm++;
    }

    cout << "Size SEQ: " << seq_nn << " " << seq_tn << " " << seq_tm << endl;
    cout << "Size K2I: " << k2i_nn << " " << k2i_tn << endl;
    cout << "Size I2P: " << i2p_tm << endl;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Time filter/stats: " << time.count() << endl;
}
