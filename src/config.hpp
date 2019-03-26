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

#ifndef __SM_CONFIG_H__
#define __SM_CONFIG_H__

#include <vector>

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
    std::vector<std::string> list_normal;
    std::vector<std::string> list_tumor;
    bool check_quality;

    std::string output_path;
    std::string output_path_count;
    std::string output_path_filter;
    std::string output_path_merge;
    std::string output_path_group;
    std::string data_path;

    std::string exec;

    double false_positive_rate;
    uint64_t all_size;
    uint64_t allowed_size;

    bool enable_cache;

    // Total number of table and cache keys.
    uint64_t table_size;
    uint64_t cache_size;

    std::string conversion_mode;
    int max_conversions;

    bool prefilter;

    int export_min;
    int export_max;

    std::string annotate_input;

    std::string index_format;
    int num_indexes;
    int max_nc_a;
    int min_tc_a;
    int max_nc_b;
    int min_tc_b;
    uint64_t leads_size;

    // Maximum number of reads per kmer while filtering; kmers with more than
    // max_filter_reads associated reads are ignored.
    int max_filter_reads;

    int window_min;
    int window_len;

    // Maximum number of reads per kmer while grouping; reads from kmers with
    // more than max_group_reads reads are dropped and marked as such in the
    // groups file.
    int max_group_reads;

    // Number of high and low priority RocksDB threads.
    int num_threads_high;
    int num_threads_low;
    uint64_t block_cache_size;
    uint64_t block_size;

    void load(const std::string &filename);

    // The following additional configuration variables are easily derived
    // from other parameters, and not read from the configuration file.
    // Available here as part of the configuration for convenience.

    int stem_len;

    // Starting position within a stem used for mapping.
    int map_pos;
};

#endif
