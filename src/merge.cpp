#include "merge.hpp"

#include <stdio.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include "db.hpp"

using std::cout;
using std::endl;
using std::string;
using namespace std::placeholders;

merge::merge(const sm_config &conf) : stage(conf)
{
    _executable["run"] = std::bind(&merge::run, this);
}

void merge::run()
{
    for (auto& kv: sm::types) {
        for (auto& set: kv.second) {
            load(kv.first, set);
        }
    }
}

void merge::load(string type, string set)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    auto it_type = sm::types.find(type);
    if (it_type == sm::types.end()) {
        cout << "Failed to merge, wrong type: " << type << endl;
        exit(1);
    }

    auto it_set = it_type->second.find(set);
    if (it_set == it_type->second.end()) {
        cout << "Failed to merge, wrong set: " << type << "/" << set << endl;
        exit(1);
    }

    rocksdb::DB* db;
    rocksdb::Options options;
    options.create_if_missing = true;
    options.IncreaseParallelism(4);
    string rdb_path = _conf.output_path + "/" + type + "-" + set + ".rdb";
    std::function<void(rocksdb::DB* db, std::string set, int i)> func;

    if (type == "seq") {
        func = std::bind(&merge::load_seq, this, _1, _2, _3);
    }

    if (type == "k2i") {
        func = std::bind(&merge::load_k2i, this, _1, _2, _3);
        options.merge_operator.reset(new IDListOperator());
    }

    if (type == "i2p") {
        func = std::bind(&merge::load_i2p, this, _1, _2, _3);
        options.merge_operator.reset(new PositionsMapOperator());
    }

    cout << "Open RocksDB: " << rdb_path << endl;
    rocksdb::Status status = rocksdb::DB::Open(options, rdb_path, &db);
    assert(status.ok());

    std::vector<std::thread> threads;
    for (int i = 0; i < _conf.num_mergers; i++)
        threads.push_back(std::thread(func, db, set, i));
    cout << "Spawned " << threads.size() << " merger threads" << endl;
    for (auto& thread: threads)
        thread.join();

    delete db;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Merge run time (" << type << "/" << set << "):"
         << time.count() << endl;
}

void merge::load_i2p(rocksdb::DB* db, string set, int i)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    string file = _conf.output_path + "/filter-" + set + "." + std::to_string(i) + ".i2p";
    cout << "Merge: " << file << endl;
    std::ifstream in(file);
    if (!in.good()) {
        cout << "Failed to open: " << file << endl;
        return;
    }

    int n = 0;
    string sid;
    sm_pos_bitmap p;

    while (in >> sid >> p.a[0] >> p.a[1] >> p.b[0] >> p.b[1]) {
        string serialized;
        encode_pos(serialized, p);
        db->Merge(rocksdb::WriteOptions(), sid, serialized);
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << i << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
}

void merge::load_k2i(rocksdb::DB* db, string set, int i)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    string file = _conf.output_path + "/filter-" + set + "." + std::to_string(i) + ".k2i";
    cout << "Merge: " << file << endl;
    std::ifstream in(file);
    if (!in.good()) {
        cout << "Failed to open: " << file << endl;
        return;
    }

    int n = 0;
    string kmer;
    int len = 0;

    while (in >> kmer >> len) {
        std::stringstream s;
        for (int i = 0; i < len; i++) {
            string sid;
            in >> sid;
            s << sid << " ";
        }
        db->Merge(rocksdb::WriteOptions(), kmer, s.str());
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << i << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
}

void merge::load_seq(rocksdb::DB* db, std::string set, int i)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    string file = _conf.output_path + "/filter-" + set + "." + std::to_string(i) + ".fq.gz";
    cout << "Merge: " << file << endl;
    gzFile in = gzopen(file.c_str(), "rb");
    if (in == NULL) {
        cout << "Failed to open: " << file << " (" << errno << ")" << endl;
        return;
    }

    int n = 0;
    int len;
    kseq_t *seq = kseq_init(in);
    rocksdb::WriteBatch batch;
    while ((len = kseq_read(seq)) >= 0) {
        string id(seq->name.s);
        string s(seq->seq.s);
        batch.Put(id, s);
        if (n % 10000 == 0) {
            db->Write(rocksdb::WriteOptions(), &batch);
            batch.Clear();
        }
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << i << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
    db->Write(rocksdb::WriteOptions(), &batch);
    kseq_destroy(seq);
    gzclose(in);
}
