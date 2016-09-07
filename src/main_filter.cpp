#include <common.hpp>
#include <main_filter.hpp>
#include <process.hpp>
#include <filter.hpp>

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

moodycamel::ConcurrentQueue<string> input_queue(1000);
boost::atomic_int input_count(0);
boost::atomic<bool> process_done(false);
int map_l1[MAP_FILE_LEN] = {0};
int map_l2[MAP_FILE_LEN] = {0};
std::vector<string> set_names = { "nn", "tm", "tn" };

sm_table* tables[NUM_STORERS];
sm_cache* caches[NUM_STORERS];
folly::ProducerConsumerQueue<sm_bulk>* queues[NUM_STORERS][MAX_LOADERS];
std::mutex filter_mutex[NUM_SETS];
std::unordered_set<string> filter_reads[NUM_SETS];
std::unordered_map<string, std::pair<std::vector<uint8_t>, std::vector<uint8_t>>> filter_i2p[NUM_SETS];
std::unordered_map<string, std::unordered_set<string>> filter_k2i[NUM_SETS];

int main(int argc, char *argv[])
{
    string input_filename;
    string map_filename;
    int pid = 0;
    int num_loaders = NUM_STORERS;
    int num_filters = NUM_STORERS;
    bool disable_stats = false;
    bool disable_filter = false;

    std::ios_base::sync_with_stdio(false);

    static const char *opts = "i:m:p:l:f:h";
    static const struct option opts_long[] = {
        { "input", required_argument, NULL, 'i' },
        { "mapping", required_argument, NULL, 'm' },
        { "pid", required_argument, NULL, 'p' },
        { "loaders", required_argument, NULL, 'l' },
        { "filters", required_argument, NULL, 'f' },
        { "disable-stats", no_argument, NULL, O_DISABLE_STATS },
        { "disable-filter", no_argument, NULL, O_DISABLE_FILTER },
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
            case O_DISABLE_STATS: disable_stats = true; break;
            case O_DISABLE_FILTER: disable_filter = true; break;
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

    if (num_loaders > MAX_LOADERS) {
        cout << "Number of loaders larger than MAX_LOADERS" << endl;
        exit(1);
    }

    // Initialize prefix to process/thread mapping.
    for (string line; getline(map_file, line);) {
        std::vector<string> columns;
        boost::split(columns, line, boost::is_any_of(" "));
        uint64_t m = 0;
        memcpy(&m, columns[0].c_str(), MAP_LEN);
        hash_5c_map(m);
        map_l1[m] = atoi(columns[1].c_str());
        map_l2[m] = atoi(columns[2].c_str());
    }

    init();

    reset_input_queue(input_file);
    sm_process(pid, num_loaders, NUM_STORERS);

    if (!disable_stats) {
        sm_stats(NUM_STORERS);
    }

    // Deallocate unnecessary memory so as to allow larger filters.
    free_tables();
    rebuild_tables();

    if (!disable_filter) {
        reset_input_queue(input_file);
        sm_filter(pid, num_filters);
    }

    return 0;
}

void display_usage()
{
    cout << "Usage: sm-filter [OPTIONS] -i INPUT_FILE -m MAP_FILE" << endl;
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

void init()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    // Initialize tables and message queues.
    for (int i = 0; i < NUM_STORERS; i++) {
        for (int j = 0; j < MAX_LOADERS; j++) {
            queues[i][j] = new folly::ProducerConsumerQueue<sm_bulk>(QMSG_LEN);
        }
    }

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Init time: " << time.count() << endl;
}

void rebuild_tables()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    std::vector<std::thread> rebuilders;
    for (int i = 0; i < NUM_STORERS; i++)
        rebuilders.push_back(std::thread(rebuild_table, i));
    cout << "Spawned " << rebuilders.size() << " rebuilder threads" << endl;

    for (auto& rebuilder: rebuilders)
        rebuilder.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Rebuild time: " << time.count() << endl;
}

void rebuild_table(int sid)
{
    tables[sid] = new sm_table();
    tables[sid]->resize(TABLE_LEN);
    string file = string("table-") + std::to_string(sid) + string(".data");
    FILE* fp = fopen(file.c_str(), "r");
    tables[sid]->unserialize(sm_table::NopointerSerializer(), fp);
}

void free_tables()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    std::vector<std::thread> deallocs;
    for (int i = 0; i < NUM_STORERS; i++)
        deallocs.push_back(std::thread(free_table, i));
    cout << "Spawned " << deallocs.size() << " dealloc threads" << endl;

    for (auto& dealloc: deallocs)
        dealloc.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Dealloc time: " << time.count() << endl;
}

void free_table(int sid)
{
    caches[sid]->clear();
    caches[sid]->resize(0);
    delete caches[sid];
    tables[sid]->clear();
    tables[sid]->resize(0);
    delete tables[sid];
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
    cout << "Filtered reads (NN+-): " << filter_reads[NN].size() << endl;
    cout << "Filtered reads (TN+-): " << filter_reads[TN].size() << endl;
    cout << "Filtered reads (TM+-): " << filter_reads[TM].size() << endl;
    cout << "Filter time: " << time.count() << endl;

#ifdef PROFILE
    ProfilerStop();
#endif

    start = std::chrono::system_clock::now();

    std::vector<std::thread> fastqs;
    for (int i = 0; i < NUM_SETS; i++)
        fastqs.push_back(std::thread(sm_write_fastq, i, pid));
    cout << "Spawned " << fastqs.size() << " FASTQ writer threads" << endl;

    std::vector<std::thread> k2is;
    for (int i = 0; i < NUM_SETS; i++)
        k2is.push_back(std::thread(sm_write_k2i, i, pid));
    cout << "Spawned " << fastqs.size() << " K2I writer threads" << endl;

    std::thread i2p = std::thread(sm_write_i2p, TM, pid);
    cout << "Spawned I2P writer thread" << endl;
    i2p.join();

    for (auto& k2i: k2is)
        k2i.join();
    for (auto& fastq: fastqs)
        fastq.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Output time: " << time.count() << endl;
}

void sm_stats(int num_storers)
{
    std::map<uint64_t, uint64_t> hist;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;

    start = std::chrono::system_clock::now();

    uint64_t subs = 0;
    uint64_t subs_unique = 0;
    uint64_t subs_cache = 0;
    for (int i = 0; i < num_storers; i++) {
        uint64_t part = 0;
        uint64_t part_unique = 0;
        for (sm_table::const_iterator it = tables[i]->begin();
             it != tables[i]->end(); ++it) {
            for (int f = 0; f < 4; f++) {
                for (int l = 0; l < 4; l++) {
                    uint16_t nc = it->second.v[f][l][NORMAL_READ];
                    uint16_t tc = it->second.v[f][l][CANCER_READ];
                    uint32_t count = nc + tc;
                    if (count == 0)
                        continue;
                    uint64_t bin = 1;
                    bin = int(log2(count));
                    hist[bin]++;
                    part = part + count;
                    part_unique++;
                }
            }
        }
        cout << KMER_LEN << "-mers (part-t-" << i << "): " << part << endl;
        cout << KMER_LEN << "-mers (part-u-" << i << "): " << part_unique << endl;
        cout << KMER_LEN << "-mers (part-c-" << i << "): " << caches[i]->size() << endl;
        subs += part;
        subs_unique += part_unique;
        subs_cache += caches[i]->size();
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
    cout << KMER_LEN << "-mers (cache):  " << subs_cache << endl;
    cout << "Iteration time:   " << time.count() << endl;
}

void sm_write_fastq(int set, int pid)
{
    std::ofstream ofs;
    ofs.open("filter-" + set_names[set] + "." + std::to_string(pid) + ".fastq");
    for (std::unordered_set<string>::const_iterator it =
         filter_reads[set].begin(); it != filter_reads[set].end(); ++it) {
        ofs << *it << endl;
    }
    ofs.close();
}

void sm_write_k2i(int set, int pid)
{
    std::ofstream ofs;
    ofs.open("filter-" + set_names[set] + "." + std::to_string(pid) + ".k2i");
    for (std::unordered_map<string, std::unordered_set<string>>::const_iterator it =
         filter_k2i[set].begin(); it != filter_k2i[set].end(); ++it) {
        if (it->second.size() > MAX_K2I_READS)
            continue;
        ofs << it->first << " " << it->second.size();
        for (std::unordered_set<string>::const_iterator sit = it->second.begin();
             sit != it->second.end(); ++sit) {
            ofs << " " << *sit;
        }
        ofs << endl;
    }
    ofs.close();
}

void sm_write_i2p(int set, int pid)
{
    std::ofstream ofs;
    ofs.open("filter-" + set_names[set] + "." + std::to_string(pid) + ".i2p");
    for (std::unordered_map<string, std::pair<std::vector<uint8_t>, std::vector<uint8_t>>>::const_iterator it =
         filter_i2p[set].begin(); it != filter_i2p[set].end(); ++it) {
        ofs << it->first << " " << it->second.first.size() << " " << it->second.second.size();
        for (std::vector<uint8_t>::const_iterator sit = it->second.first.begin();
             sit != it->second.first.end(); ++sit) {
            ofs << " " << (int) *sit;
        }
        for (std::vector<uint8_t>::const_iterator sit = it->second.second.begin();
             sit != it->second.second.end(); ++sit) {
            ofs << " " << (int) *sit;
        }
        ofs << endl;
    }
    ofs.close();
}
