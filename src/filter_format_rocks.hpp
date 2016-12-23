#ifndef __SM_FILTER_FORMAT_ROCKS_H__
#define __SM_FILTER_FORMAT_ROCKS_H__

#include <string>

#include <rocksdb/db.h>

#include "common.hpp"
#include "filter_format.hpp"

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
