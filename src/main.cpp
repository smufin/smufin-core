#include "main.hpp"

#include <getopt.h>

#include <chrono>
#include <fstream>
#include <string>
#include <thread>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

#include "count.hpp"
#include "filter.hpp"

using std::cout;
using std::endl;
using std::string;

int map_l1[MAP_FILE_LEN] = {0};
int map_l2[MAP_FILE_LEN] = {0};

int main(int argc, char *argv[])
{
    sm_config conf = sm_config();
    conf.exec = "count:run";

    static const char *opts = "i:m:p:l:f:x:h";
    static const struct option opts_long[] = {
        { "input", required_argument, NULL, 'i' },
        { "mapping", required_argument, NULL, 'm' },
        { "pid", required_argument, NULL, 'p' },
        { "loaders", required_argument, NULL, 'l' },
        { "filters", required_argument, NULL, 'f' },
        { "exec", required_argument, NULL, 'x' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 },
    };

    int opt = 0;
    int opt_index;
    while (opt != -1) {
        opt = getopt_long(argc, argv, opts, opts_long, &opt_index);
        switch (opt) {
            case 'i': conf.input_file = string(optarg); break;
            case 'm': conf.map_file = string(optarg); break;
            case 'p': conf.pid = atoi(optarg); break;
            case 'l': conf.num_loaders = atoi(optarg); break;
            case 'f': conf.num_filters = atoi(optarg); break;
            case 'x': conf.exec = string(optarg); break;
            case '?': display_usage(); return 1;
            case ':': display_usage(); return 1;
            case 'h': display_usage(); return 0;
        }
    }

    std::ifstream map_stream(conf.map_file);
    if (!map_stream.good()) {
        display_usage();
        exit(1);
    }

    if (conf.num_loaders > MAX_LOADERS) {
        cout << "Number of loaders larger than MAX_LOADERS" << endl;
        exit(1);
    }

    // Initialize 5-mer prefix to partition/storer mapping.
    for (string line; getline(map_stream, line);) {
        std::vector<string> columns;
        boost::split(columns, line, boost::is_any_of(" "));
        uint64_t m = 0;
        memcpy(&m, columns[0].c_str(), MAP_LEN);
        hash_5c_map(m);
        map_l1[m] = atoi(columns[1].c_str());
        map_l2[m] = atoi(columns[2].c_str());
    }

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

    std::unordered_map<string, stage*(*)(const sm_config &conf)> registry;
    registry["count"] = &create_stage<count>;
    registry["filter"] = &create_stage<filter>;

    std::vector<std::pair<string, stage*>> pipeline;
    for (auto& name: order) {
        if (registry.find(name) != registry.end()) {
            cout << "Initialize: " << name << endl;
            stage* s = registry[name](conf);
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

    return 0;
}

void display_usage()
{
    cout << "Usage: sm [OPTIONS] -i INPUT_FILE -m MAP_FILE -x COMMANDS" << endl;
    cout << "Options:" << endl;
    cout << " -i, --input INPUT_FILE" << endl;
    cout << " -m, --mapping MAP_FILE" << endl;
    cout << " -p, --pid ID" << endl;
    cout << " -l, --loaders NUM_LOADER_THREADS" << endl;
    cout << " -f, --filters NUM_FILTER_THREADS" << endl;
    cout << " -x, --exec COMMANDS" << endl;
    cout << " -h, --help" << endl;
}

template<typename T> stage* create_stage(const sm_config &conf) {
    return new T(conf);
}
