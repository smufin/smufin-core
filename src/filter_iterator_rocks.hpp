#ifndef __SM_FILTER_ITERATOR_ROCKS_H__
#define __SM_FILTER_ITERATOR_ROCKS_H__

#include <rocksdb/db.h>

#include "common.hpp"
#include "filter_iterator.hpp"

class seq_rocks_iterator : public filter_iterator<seq_t>
{
public:
    seq_rocks_iterator(const sm_config &conf, std::string set, int pid)
        : filter_iterator<seq_t>(conf, set, pid) {};

    bool init();
    bool next();

private:
    rocksdb::Iterator* _it;
};

#endif
