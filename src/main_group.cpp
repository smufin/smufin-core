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
l2r_table* l2r[MAX_GROUPERS];

i2r_table* i2r[2];
kmer_table* k2i[2];

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

    for (int i = 0; i < MAX_GROUPERS; i++) {
        l2p[i] = new l2p_table(100000);
        l2k[i] = new l2k_table(100000);
        l2i[i] = new l2i_table(100000);
        l2r[i] = new l2r_table(100000);
    }

    for (int i = 0; i < 2; i++) {
        i2r[i] = new i2r_table();
        k2i[i] = new kmer_table();
    }

    i2r[0]->resize(420000000);
    i2r[1]->resize(500000000);
    k2i[0]->resize(30000000);
    k2i[1]->resize(60000000);

    std::vector<string> sets = {"nn", "tn", "tm"};

    string dir;
    rocksdb::Status status;

    rocksdb::DB* i2p_db;
    rocksdb::Options options;
    options.WAL_ttl_seconds = 0;
    options.WAL_size_limit_MB = 0;
    options.merge_operator.reset(new PositionsMapOperator());
    dir = rdb + "/i2p-" + sets[TM] + ".rdb";
    cout << "Open RocksDB: " << dir << endl;
    status = rocksdb::DB::OpenForReadOnly(options, dir, &i2p_db);
    assert(status.ok());

    rocksdb::DB* i2r_tm_db;
    rocksdb::Options options2;
    options2.WAL_ttl_seconds = 0;
    options2.WAL_size_limit_MB = 0;
    dir = rdb + "/seq-" + sets[TM] + ".rdb";
    cout << "Open RocksDB: " << dir << endl;
    status = rocksdb::DB::OpenForReadOnly(options2, dir, &i2r_tm_db);
    assert(status.ok());

    string sid;
    sm_pos_bitmap p;

    int num_all = 0;
    int num_exist = 0;
    int num_map = 0;
    int num_match_a = 0;
    int num_match_b = 0;

    std::chrono::time_point<std::chrono::system_clock> istart, iend;
    std::chrono::duration<double> itime;
    istart = std::chrono::system_clock::now();

    // 1. Iterate through candidates.

    rocksdb::Iterator* it = i2p_db->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        num_all++;

        sid = it->key().ToString();
        p = decode_pos(it->value().data());

        std::vector<int> a_pos;
        std::vector<int> b_pos;

        string read;
        rocksdb::Status status;
        status = i2r_tm_db->Get(rocksdb::ReadOptions(), sid, &read);
        if (!status.ok())
            continue;

        if (num_all % 100000 == 0) {
            iend = std::chrono::system_clock::now();
            itime = iend - istart;
            cout << "S: " << num_all << " " << itime.count() << endl;
            istart = std::chrono::system_clock::now();
        }

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
            select_candidate(gid, sid, read, read, a_pos, 0);
            num_match_a++;
        }

        if (b_len >= KMIN && b_len <= KMAX && a_len > 0 && match_window(b_pos)) {
            char buf[RMAX + 1];
            copy(read.begin(), read.end(), buf);
            buf[read_length] = '\0';
            rrevcomp(buf, read_length);
            string directed_read = string(buf);
            select_candidate(gid, sid, read, directed_read, b_pos, 1);
            num_match_b++;
        }
    }

    delete it;
    delete i2p_db;
    delete i2r_tm_db;

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

    cout << "Number of kmers (0): " << k2i[0]->size() << endl;
    cout << "Number of kmers (1): " << k2i[1]->size() << endl;

    // 2. Iterate K2I retrieving all k-mers seen in candidate positions.

    for (int i = 0; i < 2; i++) {
        start = std::chrono::system_clock::now();

        rocksdb::DB* k2i_db;
        rocksdb::Options options;
        options.WAL_ttl_seconds = 0;
        options.WAL_size_limit_MB = 0;
        options.merge_operator.reset(new IDListOperator());
        dir = rdb + "/k2i-" + sets[i] + ".rdb";
        cout << "Open RocksDB: " << dir << endl;
        status = rocksdb::DB::OpenForReadOnly(options, dir, &k2i_db);
        assert(status.ok());

        int num_seen = 0;
        int num_kmer = 0;
        istart = std::chrono::system_clock::now();
        it = k2i_db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) {
            num_kmer++;

            string kmer = it->key().ToString();
            string list = it->value().ToString();
            kmer_table::const_iterator kit = k2i[i]->find(kmer);
            if (kit != k2i[i]->end()) {
                (*k2i[i])[kmer] = list;
                num_seen++;
            }

            if (num_kmer % 1000000 == 0) {
                iend = std::chrono::system_clock::now();
                itime = iend - istart;
                cout << "K: " << i << " " << num_kmer << " " << num_seen << " "
                     << itime.count() << endl;
                istart = std::chrono::system_clock::now();
            }
        }

        cout << "Number of kmers seen (" << sets[i] << "): "
             << num_seen << endl;
        cout << "Number of IDs seen (" << sets[i] << "): "
             << i2r[i]->size() << endl;


        delete it;
        delete k2i_db;

        end = std::chrono::system_clock::now();
        time = end - start;
        cout << "Kmer iteration time (" << sets[i] << "): "
             << time.count() << endl;
    }

    // 3. Iterate & collect reads.

    for (int i = 0; i < 2; i++) {
        start = std::chrono::system_clock::now();

        rocksdb::DB* i2r_db;
        dir = rdb + "/seq-" + sets[i] + ".rdb";
        cout << "Open RocksDB: " << dir << endl;
        status = rocksdb::DB::OpenForReadOnly(rocksdb::Options(), dir, &i2r_db);
        assert(status.ok());

        int num_seen = 0;
        int num_read = 0;
        istart = std::chrono::system_clock::now();
        it = i2r_db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) {
            num_read++;

            string sid = it->key().ToString();
            string seq = it->value().ToString();
            sm_read read;
            encode_read(seq, read);
            (*i2r[i])[sid] = read;

            if (num_read % 10000000 == 0) {
                iend = std::chrono::system_clock::now();
                itime = iend - istart;
                cout << "R: " << i << " " << num_read << " " << num_seen << " "
                     << itime.count() << endl;
                istart = std::chrono::system_clock::now();
            }
        }

        delete it;
        delete i2r_db;

        end = std::chrono::system_clock::now();
        time = end - start;
        cout << "Read iteration time (" << sets[i] << "): "
             << time.count() << endl;
    }

    // 4. Populate candidate groups.

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
