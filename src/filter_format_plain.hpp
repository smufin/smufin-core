#ifndef __SM_FILTER_FORMAT_PLAIN_H__
#define __SM_FILTER_FORMAT_PLAIN_H__

#include <mutex>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include "common.hpp"
#include "filter_format.hpp"

// An implementation of a filter_format that creates filtering indexes in
// memory using standard maps, and mutexes to ensure atomicity, generating and
// dumping these maps as plain space-separated text files to disk.
//
// For a particular partition P, the following 7 files are generated:
//  - filter-seq-{nn,tn,tm}.P.txt
//  - filter-k2i-{nn,tn,tm}.P.txt
//  - filter-i2p-tm.P.txt
//
// SEQ index files are periodically flushed to disk during filter execution so
// as to lower memory consumption, while K2I and I2P are only flushed at the
// end.
class filter_format_plain : public filter_format
{
public:
    filter_format_plain(const sm_config &conf) : filter_format(conf) {};

    void update(kseq_t *seq, int pos, bool rev, char kmer[], sm_idx_set set);
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
