#ifndef __SM_FILTER_H__
#define __SM_FILTER_H__

#include <string>

#include "common.hpp"
#include "count.hpp"
#include "filter_format.hpp"
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
    moodycamel::ConcurrentQueue<std::string> _input_queue;
    std::atomic<int> _input_count{0};

    const count* _count;

    filter_format* _format;

    void load(int fid);
    void load_file(int fid, std::string file);

    void filter_normal(int fid, kseq_t *seq, const char *sub, int len);
    void filter_cancer(int fid, kseq_t *seq, const char *sub, int len);

    int get_value(int fid, char kmer[], sm_table::const_iterator *it);

    inline void filter_all(int fid, kseq_t *seq, int pos, bool rev,
                           char kmer[], sm_set set);
    inline void filter_branch(int fid, kseq_t *seq, int pos, bool rev,
                              char kmer[], sm_set set);
    inline void filter_kmer(kseq_t *seq, int pos, bool rev, char kmer[],
                            uint32_t nc,uint32_t tc, uint32_t nsum,
                            uint32_t tsum, sm_set set);
};

#endif
