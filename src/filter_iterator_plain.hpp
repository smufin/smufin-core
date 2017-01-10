#ifndef __SM_FILTER_ITERATOR_PLAIN_H__
#define __SM_FILTER_ITERATOR_PLAIN_H__

#include <fstream>

#include "common.hpp"
#include "filter_iterator.hpp"

template <typename T>
class plain_iterator : public filter_iterator<T>
{
public:
    plain_iterator(const sm_config &conf, std::string set, int pid,
                   std::string type)
        : filter_iterator<T>(conf, set, pid), _type(type) {};

    bool init();

protected:
    std::ifstream _in;
    std::string _type;
};

class seq_plain_iterator : public plain_iterator<seq_t>
{
public:
    seq_plain_iterator(const sm_config &conf, std::string set, int pid)
        : plain_iterator<seq_t>(conf, set, pid, "seq") {};

    bool next();
};

class k2i_plain_iterator : public plain_iterator<k2i_t>
{
public:
    k2i_plain_iterator(const sm_config &conf, std::string set, int pid)
        : plain_iterator<k2i_t>(conf, set, pid, "k2i") {};

    bool next();
};

class i2p_plain_iterator : public plain_iterator<i2p_t>
{
public:
    i2p_plain_iterator(const sm_config &conf, std::string set, int pid)
        : plain_iterator<i2p_t>(conf, set, pid, "i2p") {};

    bool next();
};

#endif
