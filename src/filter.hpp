#ifndef __SM_FILTER_H__
#define __SM_FILTER_H__

#include <string>

#include "common.hpp"
#include "count.hpp"
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

    const char _alpha[4];
    const std::array<std::string, NUM_SETS> _sets;

    const count* _count;

    std::mutex _mutex[NUM_SETS];
    sm_ids _ids[NUM_SETS];
    sm_reads _reads[NUM_SETS];
    sm_i2p _i2p[NUM_SETS];
    sm_k2i _k2i[NUM_SETS];

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

    void write_fastq(int set);
    void write_k2i(int set);
    void write_i2p(int set);
};

#endif
