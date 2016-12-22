#ifndef __SM_FILTER_FORMAT_H__
#define __SM_FILTER_FORMAT_H__

#include "common.hpp"

class filter_format
{
public:
    filter_format(const sm_config &conf) : _conf(conf) {};

    virtual void update(kseq_t *seq, int pos, bool rev, char kmer[],
                        sm_set set) = 0;
    virtual bool flush() = 0;
    virtual void dump() = 0;
    virtual void stats() = 0;

protected:
    const sm_config &_conf;
};

#endif
