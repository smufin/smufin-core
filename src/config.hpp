#ifndef __SM_CONFIG_H__
#define __SM_CONFIG_H__

struct sm_config {
    int pid;
    int num_partitions;
    int num_loaders;
    int num_storers;
    int num_filters;
    int num_mergers;
    int num_groupers;

    std::string input_file;
    std::string output_path;
    std::string data_path;

    std::string exec;

    // Total number of table and cache keys.
    uint64_t table_size;
    uint64_t cache_size;

    int max_nc;
    int min_tc;

    // Maximum number of reads per kmer. During filtering, kmers with more
    // than max_k2i_reads associated reads are ignored.
    int max_k2i_reads;

    int window_min;
    int window_len;

    void load(const std::string &filename);
};

#endif
