#ifndef __SM_CONFIG_H__
#define __SM_CONFIG_H__

#include "common.hpp"

struct sm_config {
    int k;
    int pid;
    int num_partitions;
    int num_loaders;
    int num_storers;
    int num_filters;
    int num_mergers;
    int num_groupers;

    std::string input_format;
    std::string input_normal;
    std::string input_tumor;
    std::string output_path;
    std::string data_path;

    std::string exec;

    double false_positive_rate;
    uint64_t all_size;
    uint64_t allowed_size;

    bool enable_cache;

    // Total number of table and cache keys.
    uint64_t table_size;
    uint64_t cache_size;

    int export_min;
    int export_max;

    std::string filter_format;
    int max_nc;
    int min_tc;

    // Maximum number of reads per kmer while filtering; kmers with more than
    // max_filter_reads associated reads are ignored.
    int max_filter_reads;

    int window_min;
    int window_len;

    // Maximum number of reads per kmer while grouping; reads from kmers with
    // more than max_group_reads reads are dropped and marked as such in the
    // groups file.
    int max_group_reads;

    void load(const std::string &filename);
};

#endif
