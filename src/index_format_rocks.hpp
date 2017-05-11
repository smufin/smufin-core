#ifndef __SM_INDEX_FORMAT_ROCKS_H__
#define __SM_INDEX_FORMAT_ROCKS_H__

#include <string>

#include "common.hpp"
#include "db.hpp"
#include "index_format.hpp"

#define MAX_INDEXES 48

// Implementation of a index_format that creates filtering indexes using
// RocksDB-backed database.
//
// For a particular partition P, the following 6 RocksDB databases are
// created in the output directory:
//  - index-seq-{nn,tn,tm}.P.rdb
//  - index-k2i-{nn,tn}.P.rdb
//  - index-i2p-tm.P.rdb
//
// The flush() method is empty since RocksDB already deals with disk
// synchronization internally. On the other hand, dump() is used to force a
// compaction from L0 to L1.
class index_format_rocks : public index_format
{
public:
    index_format_rocks(const sm_config &conf);

    void update(int fid, const sm_read *read, int pos, bool rev, char kmer[],
                sm_idx_set set);
    bool flush() {};
    void dump();
    void stats();

private:
    rdb_handle _seq[NUM_SETS][MAX_INDEXES];
    rdb_handle _k2i[2][MAX_INDEXES];
    rdb_handle _i2p[MAX_INDEXES];

    void compact(rocksdb::DB* db);
};

#endif
