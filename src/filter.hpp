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
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2019
 */

#ifndef __SM_FILTER_H__
#define __SM_FILTER_H__

#include <string>

#include "common.hpp"
#include "count.hpp"
#include "index_format.hpp"
#include "input.hpp"
#include "stage.hpp"

class filter : public stage
{
public:
    filter(const sm_config &conf);
    void chain(const stage* prev);
    void run();
    void dump();
    void stats();

    // Filter condition to evaluate a kmer (A) and it reverse-complement (B).
    static inline bool condition(const sm_config &conf, uint32_t na,
                                 uint32_t ta, uint32_t nb, uint32_t tb)
    {
        if (ta >= conf.min_tc_a && na <= conf.max_nc_a &&
            tb >= conf.min_tc_b && nb <= conf.max_nc_b) {
            return true;
        }
        return false;
    }

    // Filter condition to pre-evaluate a kmer. This condition can be used
    // when prefiltering individual stems without roots: if a kmer doesn't
    // meet the condition in at least one direction, it can be guaranteed that
    // its reverse-complement won't meet the full condition either.
    static inline bool condition(const sm_config &conf, uint32_t n, uint32_t t)
    {
        if ((t >= conf.min_tc_a && n <= conf.max_nc_a) ||
            (t >= conf.min_tc_b && n <= conf.max_nc_b)) {
            return true;
        }
        return false;
    }

private:
    input_queue* _input_queue;

    const count* _count;

    index_format* _format;

    void load(int fid);
    void load_chunk(int fid, const sm_chunk &chunk);

    void filter_normal(int fid, const sm_read *read, const char *sub, int len);
    void filter_cancer(int fid, const sm_read *read, const char *sub, int len);

    int get_value(char root[], sm_root_table::const_iterator *it);

    inline void filter_all(int fid, const sm_read *read, int pos, char kmer[],
                           sm_dir dir, int order, const sm_root &counts,
                           sm_idx_set set);
    inline void filter_branch(int fid, const sm_read *read, int pos,
                              char kmer[], sm_dir dir, int order,
                              const sm_root &counts, sm_idx_set set);
    inline void filter_kmer(int fid, const sm_read *read, int pos, char kmer[],
                            sm_dir dir, uint32_t na, uint32_t ta,
                            uint32_t nb, uint32_t tb, sm_idx_set set);
};


#endif
