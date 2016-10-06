#include <getopt.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <rocksdb/db.h>
#include <rocksdb/cache.h>
#include <rocksdb/filter_policy.h>
#include <rocksdb/table.h>

#include "db.hpp"
#include "common.hpp"
#include "group.hpp"

using std::cout;
using std::endl;
using std::string;

int map_l1[MAP_FILE_LEN] = {0};
int map_l2[MAP_FILE_LEN] = {0};

l2p_table* l2p[MAX_GROUPERS];
l2k_table* l2k[MAX_GROUPERS];
l2i_table* l2i[MAX_GROUPERS];

rocksdb::DB* i2r[NUM_SETS];
rocksdb::DB* k2i[NUM_SETS];

void display_usage()
{
    cout << "Usage: sm-group [OPTIONS] -r RDB_PATH" << endl;
    cout << "Options:" << endl;
    cout << " -r, --rdb RDB_PATH" << endl;
    cout << " -m, --mapping MAP_FILE" << endl;
    cout << " -g, --groupers NUM_THREADS" << endl;
    cout << " -p, --pid PID" << endl;
    cout << " -h, --help" << endl;
}

int main(int argc, char *argv[])
{
    string rdb;
    string map_filename;
    int pid = 0;
    int num_groupers = 1;

    static const char *opts = "r:m:p:g:h";
    static const struct option opts_long[] = {
        { "rdb", required_argument, NULL, 'r' },
        { "mapping", required_argument, NULL, 'm' },
        { "pid", required_argument, NULL, 'p' },
        { "groupers", required_argument, NULL, 'g' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 },
    };

    int opt = 0;
    int opt_index;
    while (opt != -1) {
        opt = getopt_long(argc, argv, opts, opts_long, &opt_index);
        switch (opt) {
            case 'r': rdb = string(optarg); break;
            case 'm': map_filename = string(optarg); break;
            case 'p': pid = atoi(optarg); break;
            case 'g': num_groupers = atoi(optarg); break;
            case 'h':
                display_usage();
                return 0;
        }
    }

    std::ifstream map_file(map_filename);
    if (!map_file.good()) {
        display_usage();
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

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    std::vector<string> sets = {"nn", "tn", "tm"};

    string dir;
    rocksdb::Status status;

    for (int i = 0; i < NUM_SETS; i++) {
        rocksdb::Options options;
        rocksdb::BlockBasedTableOptions toptions;

        toptions.block_cache = rocksdb::NewLRUCache(8000UL * 1024 * 1024);
        toptions.filter_policy.reset(rocksdb::NewBloomFilterPolicy(10));
        options.table_factory.reset(rocksdb::NewBlockBasedTableFactory(toptions));

        dir = rdb + "/seq-" + sets[i] + ".rdb";
        cout << "Open RocksDB: " << dir << endl;
        status = rocksdb::DB::OpenForReadOnly(options, dir, &i2r[i]);
        assert(status.ok());
    }

    for (int i = 0; i < NUM_SETS; i++) {
        rocksdb::Options options;
        rocksdb::BlockBasedTableOptions toptions;

        toptions.block_cache = rocksdb::NewLRUCache(8000UL * 1024 * 1024);
        toptions.filter_policy.reset(rocksdb::NewBloomFilterPolicy(10));
        options.table_factory.reset(rocksdb::NewBlockBasedTableFactory(toptions));

        options.merge_operator.reset(new IDListOperator());

        dir = rdb + "/k2i-" + sets[i] + ".rdb";
        cout << "Open RocksDB: " << dir << endl;
        status = rocksdb::DB::OpenForReadOnly(options, dir, &k2i[i]);
        assert(status.ok());
    }

    rocksdb::DB* i2p;
    rocksdb::Options options;
    options.merge_operator.reset(new PositionsMapOperator());
    dir = rdb + "/i2p-" + sets[TM] + ".rdb";
    cout << "Open RocksDB: " << dir << endl;
    status = rocksdb::DB::OpenForReadOnly(options, dir, &i2p);
    assert(status.ok());

    for (int i = 0; i < MAX_GROUPERS; i++) {
        l2p[i] = new l2p_table();
        l2k[i] = new l2k_table();
        l2i[i] = new l2i_table();
    }

    string sid;
    sm_pos_bitmap p;

    int num_all = 0;
    int num_exist = 0;
    int num_map = 0;
    int num_match_a = 0;
    int num_match_b = 0;

    rocksdb::Iterator* it = i2p->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        num_all++;

        sid = it->key().ToString();
        p = decode_pos(it->value().data());

        std::vector<int> a_pos;
        std::vector<int> b_pos;

        string read;
        rocksdb::Status status;
        status = i2r[TM]->Get(rocksdb::ReadOptions(), sid, &read);
        if (!status.ok())
            continue;

        num_exist++;

        int read_length = read.size();

        if (read_length == 0)
            continue;

        if (read.find("N") != std::string::npos)
            continue;

        uint64_t m = 0;
        string sub = read.substr(0, MAP_LEN);
        memcpy(&m, sub.c_str(), MAP_LEN);
        hash_5c_map(m);
        if (map_l1[m] != pid)
            continue;
        int gid = map_l2[m];

        num_map++;

        get_positions_a(p.a, &a_pos);
        get_positions_b(p.b, &b_pos, read_length);
        std::reverse(b_pos.begin(), b_pos.end());

        int a_len = a_pos.size();
        int b_len = b_pos.size();

        if (a_len >= KMIN && a_len <= KMAX && b_len > 0 && match_window(a_pos)) {
            select_candidate(gid, sid, read, a_pos, 0);
            num_match_a++;
        }

        if (b_len >= KMIN && b_len <= KMAX && a_len > 0 && match_window(b_pos)) {
            char buf[RMAX + 1];
            copy(read.begin(), read.end(), buf);
            buf[read_length] = '\0';
            rrevcomp(buf, read_length);
            select_candidate(gid, sid, string(buf), b_pos, 1);
            num_match_b++;
        }
    }

    delete it;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Candidate selection time: " << time.count() << endl;

    cout << "Number of iterated I2P: " << num_all << endl;
    cout << "Number of existing I2P: " << num_exist << endl;
    cout << "Number of mapped I2P: " << num_map << endl;
    cout << "Number of A-matching I2P: " << num_match_a << endl;
    cout << "Number of B-matching I2P: " << num_match_b << endl;
    for (int i = 0; i < MAX_GROUPERS; i++) {
        cout << "Number of candidates (" << std::to_string(i) << "): "
             << l2k[i]->size() << endl;
    }

    start = std::chrono::system_clock::now();

    std::vector<std::thread> populators;
    for (int i = 0; i < num_groupers; i++)
        populators.push_back(std::thread(populate, pid, i));
    cout << "Spawned " << populators.size() << " populator threads" << endl;
    for (auto& populator: populators)
        populator.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Populate candidates time: " << time.count() << endl;
}
