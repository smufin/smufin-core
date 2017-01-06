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

#endif
