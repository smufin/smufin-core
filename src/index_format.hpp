#ifndef __SM_INDEX_FORMAT_H__
#define __SM_INDEX_FORMAT_H__

#include "common.hpp"
#include "input.hpp"

// Abstract class that represents an interface to generate filter indexes.
// Filter indexes must include mappings of sequence IDs to sequences (SEQ),
// kmers to sequence IDs (K2I), and sequence IDs to positions of candidate
// kmers (I2P). Indexes also need to categorize sequences based on the
// following sets: normal sample (NN), non-mutated tumoral sample (TN), and
// mutated tumoral sample (TM).
class index_format
{
public:
    index_format(const sm_config &conf) : _conf(conf) {};

    // Main method to add a particular position/kmer of a sequence to the
    // filter indexes.
    virtual void update(int fid, const sm_read *read, int pos, bool rev,
                        char kmer[], sm_idx_set set) = 0;
    virtual bool flush() = 0;
    virtual void dump() = 0;
    virtual void stats() = 0;

    template<typename T> static index_format* create(const sm_config &conf)
    {
        return new T(conf);
    };

protected:
    const sm_config &_conf;
};

typedef std::function<index_format*(const sm_config &conf)> index_format_s;

#endif
