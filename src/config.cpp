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

#include "config.hpp"

#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "registry.hpp"

using std::cout;
using std::endl;
using std::string;

void sm_config::load(const string &filename)
{
    boost::property_tree::ptree tree;
    std::ifstream file(filename);
    if (!file.good()) {
        cout << "Failed to load config file " << filename << endl;
        exit(1);
    }

    boost::property_tree::read_ini(filename, tree);

    k = tree.get<int>("core.k", 30);
    pid = tree.get<int>("core.pid", 0);
    num_partitions = tree.get<int>("core.num-partitions", 1);
    num_loaders = tree.get<int>("core.num-loaders", 1);
    num_storers = tree.get<int>("core.num-storers", 1);
    num_filters = tree.get<int>("core.num-filters", 1);
    num_mergers = tree.get<int>("core.num-mergers", 1);
    num_groupers = tree.get<int>("core.num-groupers", 1);

    input_format = tree.get<string>("core.input-format", "fastq");
    input_normal = tree.get<string>("core.input-normal", "");
    input_tumor = tree.get<string>("core.input-tumor", "");
    check_quality = tree.get<bool>("core.check-quality", true);

    output_path = tree.get<string>("core.output", "");
    output_path_count = tree.get<string>("count.output", output_path);
    output_path_filter = tree.get<string>("filter.output", output_path);
    output_path_merge = tree.get<string>("merge.output", output_path);
    output_path_group = tree.get<string>("group.output", output_path);

    data_path = tree.get<string>("core.data", "data");

    exec = tree.get<string>("core.exec", "count:run,stats");

    false_positive_rate = tree.get<double>("prune.false-positive-rate", 0.05);
    all_size = tree.get<uint64_t>("prune.all-size", 100000000000);
    allowed_size = tree.get<uint64_t>("prune.allowed-size", 10000000000);

    enable_cache = tree.get<bool>("count.enable-cache", true);
    table_size = tree.get<uint64_t>("count.table-size", 12800000000);
    cache_size = tree.get<uint64_t>("count.cache-size", 106240000000);
    export_min = tree.get<int>("count.export-min", 29);
    export_max = tree.get<int>("count.export-max", 31);
    annotate_input = tree.get<string>("count.annotate-input", "");
    max_conversions = tree.get<int>("count.max-conversions", num_storers);

    index_format = tree.get<string>("filter.index-format", "plain");
    num_indexes = tree.get<int>("filter.num-indexes", 1);
    max_nc_a = tree.get<int>("filter.max-normal-count-a", 1);
    min_tc_a = tree.get<int>("filter.min-tumor-count-a", 4);
    max_nc_b = tree.get<int>("filter.max-normal-count-b", 1);
    min_tc_b = tree.get<int>("filter.min-tumor-count-b", 4);
    max_filter_reads = tree.get<int>("filter.max-reads", 2000);

    window_min = tree.get<int>("group.window-min", 7);
    window_len = tree.get<int>("group.window-len", 10);
    max_group_reads = tree.get<int>("group.max-reads", 500);
    leads_size = tree.get<uint64_t>("group.leads-size", 12800000);

    num_threads_high = tree.get<int>("rocks.num-threads-high", 1);
    num_threads_low = tree.get<int>("rocks.num-threads-low", 1);
    block_cache_size = tree.get<uint64_t>("rocks.block-cache-size", 536870912);
    block_size = tree.get<uint64_t>("rocks.block-size", 4096);

    stem_len = k - 2;
    map_pos = (stem_len - MAP_LEN) / 2;

    if (sm::input_queues.find(input_format) == sm::input_queues.end()) {
        cout << "Invalid input format " << input_format << endl;
        exit(1);
    }

    if (sm::formats.find(index_format) == sm::formats.end()) {
        cout << "Invalid filter format " << index_format << endl;
        exit(1);
    }
}
