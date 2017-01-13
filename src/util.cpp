#include "util.hpp"

#include <sstream>

#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::string;

// Specialization to spawn N threads, passing a list of sequential integers:
// [0..N).
void spawn(string name, std::function<void(int)> func, int n)
{
    std::vector<int> list;
    for (int i = 0; i < n; i++)
        list.push_back(i);
    spawn<int>(name, func, list);
}

// Initialize a two-level mapping of sizes `n1' and `n2' respectively.
void init_mapping(const sm_config &conf, int n1, int n2, int l1[], int l2[])
{
    std::ostringstream map_file;
    map_file << conf.data_path << "/maps/5-" << n1 << "-" << n2;
    std::ifstream map_stream(map_file.str());
    if (!map_stream.good()) {
        cout << "Failed to load mapping " << map_file.str() << endl;
        cout << "Wrong data path and/or number of partitions/threads" << endl;
        exit(1);
    }

    // Initialize 5-mer prefix to partition/storer mapping.
    cout << "Load mapping: " << map_file.str() << endl;
    for (string line; getline(map_stream, line);) {
        std::vector<string> columns;
        boost::split(columns, line, boost::is_any_of(" "));
        uint64_t m = 0;
        memcpy(&m, columns[0].c_str(), MAP_LEN);
        hash_5c_map(m);
        l1[m] = atoi(columns[1].c_str());
        l2[m] = atoi(columns[2].c_str());
    }
}

// Estimate size in GB of a sparsehash table of n elements, with keys of size
// k and values of size v.
float estimate_sparse(uint64_t n, size_t k, size_t v)
{
    return (n * (k + v + 1)) / 1024 / 1024 / 1024;
}
