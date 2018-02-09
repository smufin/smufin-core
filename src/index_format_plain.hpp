/*
 * Copyright © 2015-2018 Barcelona Supercomputing Center (BSC)
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

#ifndef __SM_INDEX_FORMAT_PLAIN_H__
#define __SM_INDEX_FORMAT_PLAIN_H__

#include <mutex>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include "common.hpp"
#include "index_format.hpp"

// An implementation of a index_format that creates filtering indexes in
// memory using standard maps, and mutexes to ensure atomicity, generating and
// dumping these maps as plain space-separated text files to disk.
//
// For a particular partition P, the following 6 files are generated:
//  - index-seq-{nn,tn,tm}.P.txt
//  - index-k2i-{nn,tn}.P.txt
//  - index-i2p-tm.P.txt
//
// SEQ index files are periodically flushed to disk during filter execution so
// as to lower memory consumption, while K2I and I2P are only flushed at the
// end.
class index_format_plain : public index_format
{
public:
    index_format_plain(const sm_config &conf) : index_format(conf) {};

    void update(int fid, const sm_read *read, int pos, char kmer[], sm_dir dir,
                sm_idx_set set);
    bool flush();
    void dump();
    void stats();

private:
    std::mutex _mutex[NUM_SETS];
    std::unordered_set<std::string> _ids[NUM_SETS];
    std::unordered_set<std::string> _seq[NUM_SETS];
    std::unordered_map<std::string, std::unordered_set<std::string>> _k2i[2];
    std::unordered_map<std::string, sm_pos_bitmap> _i2p;

    void write_seq(sm_idx_set set);
    void write_k2i(sm_idx_set set);
    void write_i2p(sm_idx_set set);
};

#endif
