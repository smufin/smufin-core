#ifndef __SM_FILTER_ITERATOR_ROCKS_H__
#define __SM_FILTER_ITERATOR_ROCKS_H__

#include <rocksdb/db.h>

#include "common.hpp"
#include "filter_iterator.hpp"

template <typename T>
class rocks_iterator : public filter_iterator<T>
{
public:
    rocks_iterator(const sm_config &conf, sm_idx_set set, int pid,
                   sm_idx_type type)
        : filter_iterator<T>(conf, set, pid), _type(type) {};

    bool init();

protected:
    rocksdb::Iterator* _it;
    sm_idx_type _type;
};

class seq_rocks_iterator : public rocks_iterator<seq_t>
{
public:
    seq_rocks_iterator(const sm_config &conf, sm_idx_set set, int pid)
        : rocks_iterator<seq_t>(conf, set, pid, SEQ) {};
    bool next();
};

class k2i_rocks_iterator : public rocks_iterator<k2i_t>
{
public:
    k2i_rocks_iterator(const sm_config &conf, sm_idx_set set, int pid)
        : rocks_iterator<k2i_t>(conf, set, pid, K2I) {};
    bool next();
};

class i2p_rocks_iterator : public rocks_iterator<i2p_t>
{
public:
    i2p_rocks_iterator(const sm_config &conf, sm_idx_set set, int pid)
        : rocks_iterator<i2p_t>(conf, set, pid, I2P) {};
    bool next();
};

#endif
