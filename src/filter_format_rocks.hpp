#ifndef __SM_FILTER_FORMAT_ROCKS_H__
#define __SM_FILTER_FORMAT_ROCKS_H__

#include <string>

#include <rocksdb/db.h>

#include "common.hpp"
#include "filter_format.hpp"

// Implementation of a filter_format that creates filtering indexes using
// RocksDB-backed database.
//
// For a particular partition P, the following 7 RocksDB databases are
// created in the output directory:
//  - filter-seq-{nn,tn,tm}.P.rdb
//  - filter-k2i-{nn,tn,tm}.P.rdb
//  - filter-i2p-tm.P.rdb
//
// The flush() and dump() methods are empty since RocksDB already deals with
// disk synchronization internally.
class filter_format_rocks : public filter_format
{
public:
    filter_format_rocks(const sm_config &conf);

    void update(kseq_t *seq, int pos, bool rev, char kmer[], sm_set set);
    bool flush() {};
    void dump() {};
    void stats();

private:
    rocksdb::DB* _i2p;
    rocksdb::DB* _seq[NUM_SETS];
    rocksdb::DB* _k2i[NUM_SETS];
};

#endif
