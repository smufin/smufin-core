#ifndef __SM_FILTER_ITERATOR_ROCKS_H__
#define __SM_FILTER_ITERATOR_ROCKS_H__

#include <rocksdb/db.h>

#include "common.hpp"
#include "filter_iterator.hpp"

template <typename T>
class rocks_iterator : public filter_iterator<T>
{
public:
    rocks_iterator(const sm_config &conf, std::string set, int pid,
                   std::string type)
        : filter_iterator<T>(conf, set, pid), _type(type) {};

    bool init();

protected:
    rocksdb::Iterator* _it;
    std::string _type;
};

class seq_rocks_iterator : public rocks_iterator<seq_t>
{
public:
    seq_rocks_iterator(const sm_config &conf, std::string set, int pid)
        : rocks_iterator<seq_t>(conf, set, pid, "seq") {};
    bool next();
};

class k2i_rocks_iterator : public rocks_iterator<k2i_t>
{
public:
    k2i_rocks_iterator(const sm_config &conf, std::string set, int pid)
        : rocks_iterator<k2i_t>(conf, set, pid, "k2i") {};
    bool next();
};

class i2p_rocks_iterator : public rocks_iterator<i2p_t>
{
public:
    i2p_rocks_iterator(const sm_config &conf, std::string set, int pid)
        : rocks_iterator<i2p_t>(conf, set, pid, "i2p") {};
    bool next();
};

#endif
