#ifndef __SM_FILTER_FORMAT_PLAIN_H__
#define __SM_FILTER_FORMAT_PLAIN_H__

#include <mutex>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include "common.hpp"
#include "filter_format.hpp"

class filter_format_plain : public filter_format
{
public:
    filter_format_plain(const sm_config &conf) : filter_format(conf) {};

    void update(kseq_t *seq, int pos, bool rev, char kmer[], sm_set set);
    bool flush();
    void dump();
    void stats();

private:
    std::mutex _mutex[NUM_SETS];
    std::unordered_set<std::string> _ids[NUM_SETS];
    std::unordered_set<std::string> _seq[NUM_SETS];
    std::unordered_map<std::string, sm_pos_bitmap> _i2p[NUM_SETS];
    std::unordered_map<std::string, std::unordered_set<std::string>>
        _k2i[NUM_SETS];

    void write_seq(int set);
    void write_k2i(int set);
    void write_i2p(int set);
};

#endif
