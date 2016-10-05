#include <getopt.h>

#include <chrono>
#include <iostream>
#include <set>
#include <thread>
#include <vector>

#include <rocksdb/db.h>

#include "db.hpp"
#include "merge.hpp"

using std::cout;
using std::endl;
using std::string;

int main(int argc, char *argv[])
{
    string input;
    string output;
    string type;
    string set;
    int num_threads = 1;

    static const char *opts = "i:o:t:s:n:";
    static const struct option opts_long[] = {
        { "input", required_argument, NULL, 'i' },
        { "output", required_argument, NULL, 'o' },
        { "type", required_argument, NULL, 't' },
        { "set", required_argument, NULL, 's' },
        { "nthreads", required_argument, NULL, 'n' },
        { NULL, no_argument, NULL, 0 },
    };

    int opt = 0;
    int opt_index;
    while (opt != -1) {
        opt = getopt_long(argc, argv, opts, opts_long, &opt_index);
        switch (opt) {
            case 'i': input = string(optarg); break;
            case 'o': output = string(optarg); break;
            case 't': type = string(optarg); break;
            case 's': set = string(optarg); break;
            case 'n': num_threads = atoi(optarg); break;
            case '?': return 1;
            case ':': return 1;
        }
    }

    rocksdb::DB* db;
    rocksdb::Options options;
    options.create_if_missing = true;
    options.IncreaseParallelism(4);
    string rdb_path;
    std::function<void(rocksdb::DB* db, string path, string set, int i)> func;

    std::set<string> types = {"seq", "i2p", "k2i"};
    std::set<string> sets = {"nn", "tn", "tm"};

    if (types.find(type) == types.end() || sets.find(set) == sets.end()) {
        cout << "Failed to initialize, wrong type or set" << endl;
        exit(1);
    }

    if (type == "seq") {
        rdb_path = output + "seq-" + set + ".rdb";
        func = load_seq;
    }

    if (type == "i2p") {
        rdb_path = output + "i2p-" + set + ".rdb";
        func = load_i2p;
        options.merge_operator.reset(new PositionsMapOperator());
    }

    if (type == "k2i") {
        rdb_path = output + "k2i-" + set + ".rdb";
        func = load_k2i;
        options.merge_operator.reset(new IDListOperator());
    }

    cout << "Open RocksDB: " << rdb_path << endl;
    rocksdb::Status status = rocksdb::DB::Open(options, rdb_path, &db);
    assert(status.ok());

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++)
        threads.push_back(std::thread(func, db, input, set, i));
    cout << "Spawned " << threads.size() << " threads" << endl;
    for (auto& thread: threads)
        thread.join();

    delete db;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Time: " << time.count() << endl;

    return 0;
}
