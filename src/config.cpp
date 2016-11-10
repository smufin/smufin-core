#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "common.hpp"
#include "config.hpp"

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

    pid = tree.get<int>("core.pid", 0);
    num_partitions = tree.get<int>("core.num-partitions", 1);
    num_loaders = tree.get<int>("core.num-loaders", 1);
    num_storers = tree.get<int>("core.num-storers", 1);
    num_filters = tree.get<int>("core.num-filters", 1);
    num_mergers = tree.get<int>("core.num-mergers", 1);
    num_groupers = tree.get<int>("core.num-groupers", 1);

    input_file = tree.get<string>("core.input", "");
    output_path = tree.get<string>("core.output", "");
    data_path = tree.get<string>("core.data", "data");

    exec = tree.get<string>("core.exec", "count:run,stats");

    table_size = tree.get<uint64_t>("count.table-size", 12800000000);
    cache_size = tree.get<uint64_t>("count.cache-size", 106240000000);

    max_nc = tree.get<int>("filter.max-normal-count", 1);
    min_tc = tree.get<int>("filter.min-tumor-count", 1);
    max_k2i_reads = tree.get<int>("filter.max-k2i-reads", 2000);
}
