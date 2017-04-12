#ifndef __SM_FILTER_FORMAT_ROCKS_H__
#define __SM_FILTER_FORMAT_ROCKS_H__

#include <string>

#include <rocksdb/db.h>

#include "common.hpp"
#include "filter_format.hpp"

#define MAX_FILTERS 48

// Implementation of a filter_format that creates filtering indexes using
// RocksDB-backed database.
//
// For a particular partition P, the following 6 RocksDB databases are
// created in the output directory:
//  - filter-seq-{nn,tn,tm}.P.rdb
//  - filter-k2i-{nn,tn}.P.rdb
//  - filter-i2p-tm.P.rdb
//
// The flush() method is empty since RocksDB already deals with disk
// synchronization internally. On the other hand, dump() is used to force a
// compaction from L0 to L1.
class filter_format_rocks : public filter_format
{
public:
    filter_format_rocks(const sm_config &conf);

    void update(int fid, const sm_read *read, int pos, bool rev, char kmer[],
                sm_idx_set set);
    bool flush() {};
    void dump();
    void stats();

private:
    rocksdb::DB* _seq[NUM_SETS][MAX_FILTERS];
    rocksdb::DB* _k2i[2][MAX_FILTERS];
    rocksdb::DB* _i2p[MAX_FILTERS];

    void compact(rocksdb::DB* db);
};

#endif
