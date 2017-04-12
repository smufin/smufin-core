#ifndef __SM_FILTER_ITERATOR_PLAIN_H__
#define __SM_FILTER_ITERATOR_PLAIN_H__

#include <fstream>

#include "common.hpp"
#include "filter_iterator.hpp"

template <typename T>
class plain_iterator : public filter_iterator<T>
{
public:
    plain_iterator(const sm_config &conf, sm_idx_set set, int pid, int fid,
                   sm_idx_type type)
        : filter_iterator<T>(conf, set, pid, fid), _type(type) {};

    bool init();

protected:
    std::ifstream _in;
    sm_idx_type _type;
};

class seq_plain_iterator : public plain_iterator<seq_t>
{
public:
    seq_plain_iterator(const sm_config &conf, sm_idx_set set, int pid, int fid)
        : plain_iterator<seq_t>(conf, set, pid, fid, SEQ) {};
    bool next();
};

class k2i_plain_iterator : public plain_iterator<k2i_t>
{
public:
    k2i_plain_iterator(const sm_config &conf, sm_idx_set set, int pid, int fid)
        : plain_iterator<k2i_t>(conf, set, pid, fid, K2I) {};
    bool next();
};

class i2p_plain_iterator : public plain_iterator<i2p_t>
{
public:
    i2p_plain_iterator(const sm_config &conf, sm_idx_set set, int pid, int fid)
        : plain_iterator<i2p_t>(conf, set, pid, fid, I2P) {};
    bool next();
};

#endif
