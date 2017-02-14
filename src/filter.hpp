#ifndef __SM_FILTER_H__
#define __SM_FILTER_H__

#include <string>

#include "common.hpp"
#include "count.hpp"
#include "filter_format.hpp"
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

private:
    input_queue* _input_queue;

    const count* _count;

    filter_format* _format;

    void load(int fid);
    void load_chunk(int fid, const sm_chunk &chunk);

    void filter_normal(int fid, const sm_read *read, const char *sub, int len);
    void filter_cancer(int fid, const sm_read *read, const char *sub, int len);

    int get_value(int fid, char kmer[], sm_table::const_iterator *it);

    inline void filter_all(int fid, const sm_read *read, int pos, bool rev,
                           char kmer[], sm_idx_set set);
    inline void filter_branch(int fid, const sm_read *read, int pos, bool rev,
                              char kmer[], sm_idx_set set);
    inline void filter_kmer(const sm_read *read, int pos, bool rev,
                            char kmer[], uint32_t nc, uint32_t tc,
                            uint32_t nsum, uint32_t tsum, sm_idx_set set);
};

#endif
