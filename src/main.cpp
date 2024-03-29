/*
 * Copyright © 2015-2019 Barcelona Supercomputing Center (BSC)
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

#include "main.hpp"

#include <getopt.h>

#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

#include "registry.hpp"
#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

int map_l1[MAP_FILE_LEN] = {0};
int map_l2[MAP_FILE_LEN] = {0};

int main(int argc, char *argv[])
{
    sm_config conf = sm_config();

    static const char *opts = "c:p:l:s:f:m:g:o:x:vh";
    static const struct option opts_long[] = {
        { "config", required_argument, NULL, 'c' },
        { "pid", required_argument, NULL, 'P' },
        { "partitions", required_argument, NULL, 'p' },
        { "loaders", required_argument, NULL, 'l' },
        { "storers", required_argument, NULL, 's' },
        { "filters", required_argument, NULL, 'f' },
        { "mergers", required_argument, NULL, 'm' },
        { "groupers", required_argument, NULL, 'g' },
        { "input-normal", required_argument, NULL, 'N' },
        { "input-tumor", required_argument, NULL, 'T' },
        { "output", required_argument, NULL, 'o' },
        { "exec", required_argument, NULL, 'x' },
        { "version", no_argument, NULL, 'v' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 },
    };

    int opt = 0;
    int opt_index;
    while (opt != -1) {
        opt = getopt_long(argc, argv, opts, opts_long, &opt_index);
        switch (opt) {
            case 'c': conf.load(string(optarg)); break;
            case 'P': conf.pid = atoi(optarg); break;
            case 'p': conf.num_partitions = atoi(optarg); break;
            case 'l': conf.num_loaders = atoi(optarg); break;
            case 's': conf.num_storers = atoi(optarg); break;
            case 'f': conf.num_filters = atoi(optarg); break;
            case 'm': conf.num_mergers = atoi(optarg); break;
            case 'g': conf.num_groupers = atoi(optarg); break;
            case 'N': conf.input_normal = string(optarg); break;
            case 'T': conf.input_tumor= string(optarg); break;
            case 'o':
                // Override all output paths in the configuration.
                conf.output_path = string(optarg);
                conf.output_path_count = string(optarg);
                conf.output_path_filter = string(optarg);
                conf.output_path_merge = string(optarg);
                conf.output_path_group = string(optarg);
                break;
            case 'x': conf.exec = string(optarg); break;
            case 'v': cout << VERSION << endl; return 0;
            case '?': display_usage(); return 1;
            case ':': display_usage(); return 1;
            case 'h': display_usage(); return 0;
        }
    }

    if (conf.pid >= conf.num_partitions) {
        cout << "Partition ID larger than number of partitions" << endl;
        exit(1);
    }

    conf.list_normal = expand_path(conf.input_normal);
    conf.list_tumor = expand_path(conf.input_tumor);

    cout << "Partition: " << conf.pid << " [" << conf.num_partitions
         << "]" << endl;

    init_mapping(conf, conf.num_partitions, conf.num_storers, map_l1, map_l2);

    // Parse reprogrammable execution.
    std::vector<string> commands;
    std::vector<string> order;
    std::unordered_map<string, std::vector<string>> stages;
    boost::split(commands, conf.exec, boost::is_any_of(";"));
    for (auto& command: commands) {
        size_t p = command.find_first_of(":");
        string name = command.substr(0, p);
        string list = command.substr(p + 1, string::npos);
        std::vector<string> steps;
        boost::split(steps, list, boost::is_any_of(","));
        order.push_back(name);
        stages[name] = steps;
    }

    std::vector<std::pair<string, stage*>> pipeline;
    for (auto& name: order) {
        if (sm::stages.find(name) != sm::stages.end()) {
            cout << "Initialize: " << name << endl;
            stage* s = sm::stages.at(name)(conf);
            pipeline.push_back(std::pair<string, stage*>(name, s));
        } else {
            cout << "Skip unknown stage: " << name << endl;
        }
    }

    int num_pipe = 0;
    for (auto& pipe: pipeline) {
        string name = pipe.first;
        stage* s = pipe.second;
        if (num_pipe > 0)
            pipe.second->chain(pipeline[num_pipe - 1].second);
        pipe.second->exec(name, stages[name]);
        num_pipe++;
    }

    fflush(NULL);
    _exit(0);
}

void display_usage()
{
    cout << "Usage: sm -c CONFIG [OPTIONS]" << endl;
    cout << "Options:" << endl;
    cout << " -p, --partitions NUM_PARTITIONS" << endl;
    cout << " --pid PARTITION_ID" << endl;
    cout << " -l, --loaders NUM_LOADERS" << endl;
    cout << " -s, --storers NUM_STORERS" << endl;
    cout << " -f, --filters NUM_FILTERS" << endl;
    cout << " --input-normal INPUT_FILES" << endl;
    cout << " --input-tumor INPUT_FILES" << endl;
    cout << " -o, --output OUTPUT_PATH" << endl;
    cout << " -x, --exec COMMANDS" << endl;
    cout << " -h, --help" << endl;
}
