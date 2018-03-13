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

#include "util.hpp"

#include <wordexp.h>

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
    map_file << conf.data_path << "/maps/6-" << n1 << "-" << n2;
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
        hash_6mer(m);
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

// Perform expansion of the given path, returning a vector with all matches.
std::vector<string> expand_path(string path)
{
    std::vector<string> expanded;
    wordexp_t we;
    wordexp(path.c_str(), &we, WRDE_NOCMD);
    for (int i = 0; i < we.we_wordc; i++)
        expanded.push_back(string(we.we_wordv[i]));
    wordfree(&we);
    return expanded;
}
