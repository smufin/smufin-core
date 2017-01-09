#include "merge.hpp"

#include <stdio.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>

#include "db.hpp"
#include "filter_iterator_plain.hpp"

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

    std::function<void(rocksdb::DB* db, std::string set, int i)> func;

    if (type == "seq") {
        func = std::bind(&merge::load_seq, this, _1, _2, _3);
    }

    if (type == "k2i") {
        func = std::bind(&merge::load_k2i, this, _1, _2, _3);
    }

    if (type == "i2p") {
        func = std::bind(&merge::load_i2p, this, _1, _2, _3);
    }

    rocksdb::DB* db;
    rocksdb::Options options;
    set_options_merge(options);
    set_options_type(options, type);
    string rdb = _conf.output_path + "/filter-" + type + "-" + set + ".rdb";
    cout << "Open RocksDB: " << rdb << endl;
    rocksdb::Status status = rocksdb::DB::Open(options, rdb, &db);
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
    cout << "Merge run time (" << type << "/" << set << "): "
         << time.count() << endl;
}

void merge::load_i2p(rocksdb::DB* db, string set, int pid)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    i2p_plain_iterator it(_conf, set, pid);
    if (!it.init())
        return;

    uint64_t n = 0;
    while (it.next()) {
        const i2p_t *i = it.get();
        string serialized;
        encode_pos(i->second, serialized);
        db->Merge(rocksdb::WriteOptions(), i->first, serialized);
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << pid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
}

void merge::load_k2i(rocksdb::DB* db, string set, int pid)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    k2i_plain_iterator it(_conf, set, pid);
    if (!it.init())
        return;

    uint64_t n = 0;
    while (it.next()) {
        const k2i_t *i = it.get();
        db->Merge(rocksdb::WriteOptions(), i->first, i->second);
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << pid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
}

void merge::load_seq(rocksdb::DB* db, std::string set, int pid)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    seq_plain_iterator it(_conf, set, pid);
    if (!it.init())
        return;

    uint64_t n = 0;
    rocksdb::WriteBatch batch;
    while (it.next()) {
        const seq_t *i = it.get();
        batch.Put(i->first, i->second);
        if (n % 10000 == 0) {
            db->Write(rocksdb::WriteOptions(), &batch);
            batch.Clear();
        }
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << pid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
}
