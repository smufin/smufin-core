/*
 * Copyright © 2015-2019 Barcelona Supercomputing Center (BSC)
 *
 * This file is part of SMUFIN Core. SMUFIN Core is released under the SMUFIN
 * Public License, and may not be used except in compliance with it. This file
 * is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; see the SMUFIN Public License for more details. You should have
 * received a copy of the SMUFIN Public License along with this file. If not,
 * see <https://github.com/smufin/smufin-core/blob/master/COPYING>.
 *
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2018
 */

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
    virtual void update(int fid, const sm_read *read, int pos, char kmer[],
                        sm_dir dir, sm_idx_set set) = 0;
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
