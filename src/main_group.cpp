#include <getopt.h>
#include <unistd.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <set>

#include <boost/algorithm/string.hpp>
#include <google/sparse_hash_map>
#include <rocksdb/db.h>

#include "db.hpp"
#include "common.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

#define RMAX 100
#define KMIN 0
#define KMAX 100

#define DROP 500

int map_l1[MAP_FILE_LEN] = {0};
int map_l2[MAP_FILE_LEN] = {0};

rocksdb::DB* i2r[NUM_SETS];
rocksdb::DB* k2i[NUM_SETS];

typedef std::array<std::vector<int>, 2> pos_value;
typedef std::array<std::vector<string>, 2> kmer_value;
typedef std::array<std::unordered_set<string>, 2> id_value;

typedef google::sparse_hash_map<string, pos_value> l2p_table;
typedef google::sparse_hash_map<string, kmer_value> l2k_table;
typedef google::sparse_hash_map<string, id_value> l2i_table;

typedef std::array<std::unordered_map<string, int>, 2> index_count;

l2p_table l2p;
l2k_table l2k;
l2i_table l2i;

const char comp_code[] = "ab";
const char kind_code[] = "nt";

void rrevcomp(char read[], int len)
{
    const char comp[] = "-------------------------------------------"
                        "----------------------T-G---C------N-----A";
    int c, i, j;
    for (i = 0, j = len - 1; i < j; i++, j--) {
        c = read[i];
        read[i] = comp[read[j]];
        read[j] = comp[c];
    }
}

bool match_window(std::vector<int> pos)
{
    if (pos.size() < WMIN) {
        return false;
    }

    for (int i = 0; i <= pos.size() - WMIN; i++) {
        std::vector<int> sub(pos.begin() + i, pos.begin() + i + WMIN);
        if (sub.back() - sub.front() < WLEN) {
            return true;
        }
    }

    return false;
}

void select_candidate(string sid, string seq, std::vector<int>& pos, int dir)
{
    l2p[sid][dir] = pos;
    // TODO: Switch to unordered_set instead of vector?
    std::vector<string> kmers;
    for (int p: pos) {
        kmers.push_back(seq.substr(p, KMER_LEN));
    }
    l2k[sid][dir] = kmers;
}

void populate_index(string& lid, const std::vector<string>& kmers, int kind,
                    index_count& keep, index_count& drop)
{
    for (string kmer: kmers) {
        string list;
        rocksdb::Status status;
        status = k2i[kind]->Get(rocksdb::ReadOptions(), kmer, &list);
        if (!status.ok())
            continue;

        std::unordered_set<string> sids;
        boost::split(sids, list, boost::is_any_of(" "));

        if (sids.size() > DROP) {
            drop[kind][kmer] += sids.size();
            continue;
        }

        l2i[lid][kind].insert(sids.begin(), sids.end());
        keep[kind][kmer] += sids.size();
    }
}

void get_positions_a(uint64_t bitmap[2], std::vector<int> *pos)
{
    for (int i = 0; i < 2; i++) {
        unsigned long tmp = bitmap[i];
        int offset = i * 64;
        while (tmp > 0) {
            int p = __builtin_ffsl(tmp) - 1;
            tmp &= (tmp - 1);
            pos->push_back(p + offset);
        }
    }
}

void get_positions_b(uint64_t bitmap[2], std::vector<int> *pos, int len)
{
    for (int i = 0; i < 2; i++) {
        unsigned long tmp = bitmap[i];
        int offset = i * 64;
        while (tmp > 0) {
            int p = __builtin_ffsl(tmp) - 1;
            tmp &= (tmp - 1);
            pos->push_back(len - KMER_LEN - (p + offset));
        }
    }
}

void display_usage()
{
    cout << "Usage: sm-group [OPTIONS] -i INPUT_PATH" << endl;
    cout << "Options:" << endl;
    cout << " -i, --input INPUT_PATH" << endl;
    cout << " -h, --help" << endl;
}

int main(int argc, char *argv[])
{
    string input;
    string rdb;
    string map_filename;
    int pid = 0;

    static const char *opts = "i:r:m:p:h";
    static const struct option opts_long[] = {
        { "input", required_argument, NULL, 'i' },
        { "rdb", required_argument, NULL, 'r' },
        { "mapping", required_argument, NULL, 'm' },
        { "pid", required_argument, NULL, 'p' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 },
    };

    int opt = 0;
    int opt_index;
    while (opt != -1) {
        opt = getopt_long(argc, argv, opts, opts_long, &opt_index);
        switch (opt) {
            case 'i': input = string(optarg); break;
            case 'r': rdb = string(optarg); break;
            case 'm': map_filename = string(optarg); break;
            case 'p': pid = atoi(optarg); break;
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
        dir = rdb + "/seq-" + sets[i] + ".rdb";
        cerr << "Open RocksDB: " << dir << endl;
        status = rocksdb::DB::OpenForReadOnly(options, dir, &i2r[i]);
        assert(status.ok());
    }

    for (int i = 0; i < NUM_SETS; i++) {
        rocksdb::Options options;
        options.merge_operator.reset(new IDListOperator());
        dir = rdb + "/k2i-" + sets[i] + ".rdb";
        cerr << "Open RocksDB: " << dir << endl;
        status = rocksdb::DB::OpenForReadOnly(options, dir, &k2i[i]);
        assert(status.ok());
    }

    rocksdb::DB* i2p;
    rocksdb::Options options;
    options.merge_operator.reset(new PositionsMapOperator());
    dir = rdb + "/i2p-" + sets[TM] + ".rdb";
    cerr << "Open RocksDB: " << dir << endl;
    status = rocksdb::DB::OpenForReadOnly(options, dir, &i2p);
    assert(status.ok());

    string sid;
    sm_pos_bitmap p;

    rocksdb::Iterator* it = i2p->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        sid = it->key().ToString();
        p = decode_pos(it->value().data());

        std::vector<int> a_pos;
        std::vector<int> b_pos;

        string read;
        rocksdb::Status status;
        status = i2r[TM]->Get(rocksdb::ReadOptions(), sid, &read);
        if (!status.ok())
            continue;

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

        get_positions_a(p.a, &a_pos);
        get_positions_b(p.b, &b_pos, read_length);
        std::reverse(b_pos.begin(), b_pos.end());

        int a_len = a_pos.size();
        int b_len = b_pos.size();

        if (a_len >= KMIN && a_len <= KMAX && b_len > 0 && match_window(a_pos)) {
            select_candidate(sid, read, a_pos, 0);
        }

        if (b_len >= KMIN && b_len <= KMAX && a_len > 0 && match_window(b_pos)) {
            char buf[RMAX + 1];
            copy(read.begin(), read.end(), buf);
            buf[read_length] = '\0';
            rrevcomp(buf, read_length);
            select_candidate(sid, string(buf), b_pos, 1);
        }
    }
    delete it;

    end = std::chrono::system_clock::now();
    time = end - start;
    cerr << "Candidate selection time: " << time.count() << endl;

    cout << "{";
    start = std::chrono::system_clock::now();
    bool first_group = true;
    for (l2k_table::const_iterator it = l2k.begin(); it != l2k.end(); ++it) {
        string lid = it->first;
        index_count keep;
        index_count drop;

        for (int i = 0; i < 2; i++) {
            for (string kmer: it->second[i]) {
                keep[i][kmer] = 0;
                drop[i][kmer] = 0;
            }
        }

        populate_index(lid, it->second[0], NN, keep, drop);
        populate_index(lid, it->second[0], TN, keep, drop);
        populate_index(lid, it->second[1], NN, keep, drop);
        populate_index(lid, it->second[1], TN, keep, drop);

        if (!first_group)
            cout << ",";
        first_group = false;

        string read;
        rocksdb::Status status;
        status = i2r[TM]->Get(rocksdb::ReadOptions(), lid, &read);
        if (!status.ok())
            continue;

        cout << "\"" << lid << "\":{";
        cout << "\"lead\":["
             << "\"" << lid << "\","
             << "\"" << read << "\""
             << "],";

        for (int i = 0; i < 2; i++) {
            cout << "\"pos-" << comp_code[i] << "\":[";
            bool first_pos = true;
            for (int p: l2p[lid][i]) {
                if (!first_pos)
                    cout << ",";
                first_pos = false;
                cout << p;
            }
            cout << "],";

            cout << "\"kmers-" << comp_code[i] << "\":[";
            bool first_kmer = true;
            for (string kmer: it->second[i]) {
                int kept_n = keep[0][kmer];
                int dropped_n = drop[0][kmer];
                int kept_t = keep[1][kmer];
                int dropped_t = drop[1][kmer];
                if (!first_kmer)
                    cout << ",";
                first_kmer = false;
                cout << "[\"" << kmer << "\"," << kept_n << "," << kept_t << ","
                     << dropped_n << "," << dropped_t << "]";
            }
            cout << "],";
        }

        for (int i = 0; i < 2; i++) {
            bool first_read = true;
            cout << "\"reads-" << kind_code[i] << "\":[";
            for (string sid: l2i[lid][i]) {
                if (!first_read)
                    cout << ",";
                first_read = false;
                read = "";
                status = i2r[i]->Get(rocksdb::ReadOptions(), lid, &read);
                if (!status.ok())
                    continue;
                cout << "[\"" << sid << "\",\"" << read << "\"]";
            }
            cout << ( i == 1 ? "]" : "]," );
        }

        cout << "}";
    }
    cout << "}";
    end = std::chrono::system_clock::now();
    time = end - start;
    cerr << "Populate candidates time: " << time.count() << endl;
}
