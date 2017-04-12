#ifndef __SM_MERGE_H__
#define __SM_MERGE_H__

#include <concurrentqueue.h>
#include <rocksdb/db.h>

#include "common.hpp"
#include "stage.hpp"

#include "filter_iterator.hpp"

// Signature for functions that loads index data of a particular type and set
// and a specific partition (pid) into a RocksDB database (db). E.g.
// merge::load_{seq,k2i,i2p}.
typedef std::function<void(rocksdb::DB* db, sm_idx_set set, int pid, int iid)> load_f;

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
    void load_indexes(load_f load_indexes, rocksdb::DB* db, sm_idx_set set);

    void load_seq(rocksdb::DB* db, sm_idx_set set, int pid, int iid);
    void load_k2i(rocksdb::DB* db, sm_idx_set set, int pid, int iid);
    void load_i2p(rocksdb::DB* db, sm_idx_set set, int pid, int iid);
};

#endif
