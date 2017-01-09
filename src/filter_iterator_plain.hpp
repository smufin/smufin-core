#ifndef __SM_FILTER_ITERATOR_PLAIN_H__
#define __SM_FILTER_ITERATOR_PLAIN_H__

#include <fstream>

#include "common.hpp"
#include "filter_iterator.hpp"

class seq_plain_iterator : public filter_iterator<seq_t>
{
public:
    seq_plain_iterator(const sm_config &conf, std::string set, int pid)
        : filter_iterator<seq_t>(conf, set, pid) {};

    bool init();
    bool next();

private:
    std::ifstream _in;
};

class k2i_plain_iterator : public filter_iterator<k2i_t>
{
public:
    k2i_plain_iterator(const sm_config &conf, std::string set, int pid)
        : filter_iterator<k2i_t>(conf, set, pid) {};

    bool init();
    bool next();

private:
    std::ifstream _in;
};

class i2p_plain_iterator : public filter_iterator<i2p_t>
{
public:
    i2p_plain_iterator(const sm_config &conf, std::string set, int pid)
        : filter_iterator<i2p_t>(conf, set, pid) {};

    bool init();
    bool next();

private:
    std::ifstream _in;
};

#endif
