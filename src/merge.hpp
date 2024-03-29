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

#ifndef __SM_MERGE_H__
#define __SM_MERGE_H__

#include <concurrentqueue.h>
#include <rocksdb/db.h>

#include "common.hpp"
#include "stage.hpp"

#include "db.hpp"
#include "index_iterator.hpp"

// Signature for functions that loads index data of a particular type and set
// and a specific partition (pid) into a RocksDB database (db). E.g.
// merge::load_{seq,k2i,i2p}.
typedef std::function<void(rdb_handle &rdb, sm_idx_set set, int pid,
        int iid)> load_f;

// The merge stage combines partial filtering results from multiple partitions
// into a single set of filter indexes.
class merge : public stage
{
public:
    merge(const sm_config &conf);
    void run();
    void stats();

private:
    moodycamel::ConcurrentQueue<std::pair<int, int>> _index_queue;

    void load(sm_idx_type type, sm_idx_set set);
    void load_indexes(load_f load_indexes, rdb_handle &rdb, sm_idx_set set);

    void load_seq(rdb_handle &rdb, sm_idx_set set, int pid, int iid);
    void load_k2i(rdb_handle &rdb, sm_idx_set set, int pid, int iid);
    void load_i2p(rdb_handle &rdb, sm_idx_set set, int pid, int iid);

    void to_fastq();
    void to_fastq_set(sm_idx_set set);
};

#endif
