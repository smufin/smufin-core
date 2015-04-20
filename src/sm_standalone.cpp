#include <sm_common.hpp>
#include <sm_standalone.hpp>
#include <sm_process.hpp>
#include <sm_filter.hpp>

#include <getopt.h>
#include <string>
#include <iostream>
#include <map>
#include <chrono>
#include <thread>
#include <boost/algorithm/string.hpp>
#ifdef PROFILE
#include <gperftools/profiler.h>
#endif

using std::cout;
using std::endl;
using std::string;

moodycamel::ConcurrentQueue<string> input_queue(200);
boost::atomic_int input_count(0);
boost::atomic<bool> process_done(false);
int map_l1[MAP_FILE_LEN] = {0};
int map_l2[MAP_FILE_LEN] = {0};
sm_table tables[NUM_STORERS];
folly::ProducerConsumerQueue<sm_bulk>* queues[NUM_STORERS][MAX_NUM_LOADERS];
std::unordered_set<string> filter_reads;
std::mutex filter_mutex;

int main(int argc, char *argv[])
{
    string input_filename;
    string map_filename;
    int pid = 0;
    int num_loaders = NUM_STORERS;
    int num_filters = NUM_STORERS;
    bool disable_filter = false;
    bool disable_stats = false;

    std::ios_base::sync_with_stdio(false);

    static const char *opts = "i:m:p:l:f:h";
    static const struct option opts_long[] = {
        { "input", required_argument, NULL, 'i' },
        { "mapping", required_argument, NULL, 'm' },
        { "pid", required_argument, NULL, 'p' },
        { "loaders", required_argument, NULL, 'l' },
        { "filters", required_argument, NULL, 'f' },
        { "disable-filter", no_argument, NULL, O_DISABLE_FILTER },
        { "disable-stats", no_argument, NULL, O_DISABLE_STATS },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 },
    };

    int opt = 0;
    int opt_index;
    while (opt != -1) {
        opt = getopt_long(argc, argv, opts, opts_long, &opt_index);
        switch (opt) {
            case 'i': input_filename = string(optarg); break;
            case 'm': map_filename = string(optarg); break;
            case 'p': pid = atoi(optarg); break;
            case 'l': num_loaders = atoi(optarg); break;
            case 'f': num_filters = atoi(optarg); break;
            case O_DISABLE_FILTER: disable_filter = true; break;
            case O_DISABLE_STATS: disable_stats = true; break;
            case 'h':
                display_usage();
                return 0;
        }
    }

    std::ifstream input_file(input_filename);
    std::ifstream map_file(map_filename);
    if (!input_file.good() || !map_file.good()) {
        display_usage();
        exit(1);
    }

    if (num_loaders > MAX_NUM_LOADERS) {
        cout << "Number of loaders larger than MAX_NUM_LOADERS" << endl;
        exit(1);
    }

    // Initialize prefix to process/thread mapping.
    for (string line; getline(map_file, line);) {
        std::vector<string> columns;
        boost::split(columns, line, boost::is_any_of(" "));
        uint32_t m = 0;
        memcpy(&m, columns[0].c_str(), 4);
        hash_4c_map(m);
        map_l1[m] = atoi(columns[1].c_str());
        map_l2[m] = atoi(columns[2].c_str());
    }

    // Initialize tables and message queues.
    for (int i = 0; i < NUM_STORERS; i++) {
        tables[i].resize(TABLE_LEN);
        for (int j = 0; j < MAX_NUM_LOADERS; j++) {
            queues[i][j] = new folly::ProducerConsumerQueue<sm_bulk>(QMSG_LEN);
        }
    }

    reset_input_queue(input_file);
    sm_process(pid, num_loaders, NUM_STORERS);

    if (!disable_filter) {
        reset_input_queue(input_file);
        sm_filter(pid, num_filters);
    }

    if (!disable_stats) {
        sm_stats(NUM_STORERS);
    }

    return 0;
}

void display_usage()
{
    cout << "Usage: sm-standalone [OPTIONS] -i INPUT_FILE -m MAP_FILE" << endl;
    cout << "Options:" << endl;
    cout << " -i, --input INPUT_FILE" << endl;
    cout << " -m, --mapping MAP_FILE" << endl;
    cout << " -p, --pid ID" << endl;
    cout << " -l, --loaders NUM_LOADER_THREADS" << endl;
    cout << " -f, --filters NUM_FILTER_THREADS" << endl;
    cout << " --disable-filter" << endl;
    cout << " --disable-stats" << endl;
    cout << " -h, --help" << endl;
}

void reset_input_queue(std::ifstream &input_file)
{
    input_file.clear();
    input_file.seekg(0, std::ios::beg);
    input_count = 0;
    for (string line; std::getline(input_file, line);) {
        input_queue.enqueue(line);
        input_count++;
    }
}

void sm_process(int pid, int num_loaders, int num_storers)
{
#ifdef PROFILE
    ProfilerStart("sm-process.prof");
#endif

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;

    start = std::chrono::system_clock::now();

    std::vector<std::thread> loaders;
    for (int i = 0; i < num_loaders; i++)
        loaders.push_back(std::thread(process_load, pid, i));
    cout << "Spawned " << loaders.size() << " loader threads" << endl;

    std::vector<std::thread> storers;
    for (int i = 0; i < num_storers; i++)
        storers.push_back(std::thread(process_incr, i, num_loaders));
    cout << "Spawned " << storers.size() << " storer threads" << endl;

    for (auto& loader: loaders)
        loader.join();
    process_done = true;
    for (auto& storer: storers)
        storer.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Process time: " << time.count() << endl;

#ifdef PROFILE
    ProfilerStop();
#endif
}

void sm_filter(int pid, int num_filters)
{
#ifdef PROFILE
    ProfilerStart("sm-filter.prof");
#endif

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;

    start = std::chrono::system_clock::now();

    std::vector<std::thread> filters;
    for (int i = 0; i < num_filters; i++)
        filters.push_back(std::thread(filter, pid, i));
    cout << "Spawned " << filters.size() << " filter threads" << endl;

    for (auto& filter: filters)
        filter.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Filtered reads: " << filter_reads.size() << endl;
    cout << "Filter time: " << time.count() << endl;

#ifdef PROFILE
    ProfilerStop();
#endif

    std::ofstream ofs;
    ofs.open("filtered.fastq");
    for (std::unordered_set<string>::const_iterator it = filter_reads.begin();
         it != filter_reads.end(); ++it) {
        ofs << *it << endl;
    }
}

void sm_stats(int num_storers)
{
    std::map<uint64_t, uint64_t> hist;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;

    start = std::chrono::system_clock::now();

    uint64_t subs = 0;
    uint64_t subs_unique = 0;
    for (int i = 0; i < num_storers; i++) {
        uint64_t part = 0;
        uint64_t part_unique = 0;
        for (sm_table::const_iterator it = tables[i].begin();
             it != tables[i].end(); ++it) {
            uint64_t count = it->second.first + it->second.second;
            uint64_t bin = 1;
            bin = int(log2(count));
            hist[bin]++;
            part = part + count;
            part_unique++;
        }
        cout << KMER_LEN << "-mers (part-t-" << i << "): " << part << endl;
        cout << KMER_LEN << "-mers (part-u-" << i << "): " << part_unique << endl;
        subs += part;
        subs_unique += part_unique;
    }

    end = std::chrono::system_clock::now();

    for (std::map<uint64_t, uint64_t>::const_iterator it = hist.begin();
         it != hist.end(); ++it) {
        cout << "Histo " << it->first << " " << exp2(it->first)
             << " " << it->second << endl;
    }

    time = end - start;
    cout << KMER_LEN << "-mers (total):  " << subs << endl;
    cout << KMER_LEN << "-mers (unique): " << subs_unique << endl;
    cout << "Iteration time:   " << time.count() << endl;
}
