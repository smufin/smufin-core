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

    _executable["run_seq_nn"] = std::bind(&merge::load, this, SEQ, NN);
    _executable["run_seq_tn"] = std::bind(&merge::load, this, SEQ, TN);
    _executable["run_seq_tm"] = std::bind(&merge::load, this, SEQ, TM);
    _executable["run_k2i_nn"] = std::bind(&merge::load, this, K2I, NN);
    _executable["run_k2i_tn"] = std::bind(&merge::load, this, K2I, TN);
    _executable["run_k2i_tm"] = std::bind(&merge::load, this, K2I, TM);
    _executable["run_i2p_tm"] = std::bind(&merge::load, this, I2P, TM);

    _executable["stats"] = std::bind(&merge::stats, this);
}

void merge::run()
{
    for (auto& kv: sm::indexes) {
        for (auto& set: kv.second) {
            load(kv.first, set);
        }
    }
}

// Create a RocksDB instance for a merged index, and load data from all
// partitions for a given type and set.
void merge::load(sm_idx_type type, sm_idx_set set)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    for (int i = 0; i < _conf.num_partitions; i++) {
        _partition_queue.enqueue(i);
    }

    std::map<sm_idx_type, load_f> load_map;
    load_map[SEQ] = std::bind(&merge::load_seq, this, _1, _2, _3);
    load_map[K2I] = std::bind(&merge::load_k2i, this, _1, _2, _3);
    load_map[I2P] = std::bind(&merge::load_i2p, this, _1, _2, _3);

    rocksdb::DB* db;
    open_merge(&db, _conf, type, set);

    std::vector<std::thread> threads;
    for (int i = 0; i < _conf.num_mergers; i++)
        threads.push_back(std::thread(std::bind(&merge::load_partitions, this,
                          load_map[type], db, set)));
    cout << "Spawned " << threads.size() << " merger threads" << endl;
    for (auto& thread: threads)
        thread.join();

    delete db;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Merge run time (" << sm::types[type] << "/" << sm::sets[set]
         << "): " << time.count() << endl;
}

// Load all partitions of a certain type using `load_part' and the given set
// to a RocksDB database.
void merge::load_partitions(load_f load_part, rocksdb::DB* db, sm_idx_set set)
{
    int pid;
    while (_partition_queue.try_dequeue(pid)) {
        cout << "Dequeued: " << pid << endl;
        load_part(db, set, pid);
    }
}

// Load I2P index data for a given set and partition `pid' to the database.
void merge::load_i2p(rocksdb::DB* db, sm_idx_set set, int pid)
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

// Load K2I index data for a given set and partition `pid' to the database.
void merge::load_k2i(rocksdb::DB* db, sm_idx_set set, int pid)
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

// Load SEQ index data for a given set and partition `pid' to the database.
void merge::load_seq(rocksdb::DB* db, sm_idx_set set, int pid)
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

    db->Write(rocksdb::WriteOptions(), &batch);
}

void merge::stats()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    rocksdb::DB* i2p;
    rocksdb::DB* seq[NUM_SETS];
    rocksdb::DB* k2i[NUM_SETS];

    open_merge(&i2p, _conf, I2P, TM, true);
    open_merge(&seq[NN], _conf, SEQ, NN, true);
    open_merge(&seq[TN], _conf, SEQ, TN, true);
    open_merge(&seq[TM], _conf, SEQ, TM, true);
    open_merge(&k2i[NN], _conf, K2I, NN, true);
    open_merge(&k2i[TN], _conf, K2I, TN, true);
    open_merge(&k2i[TM], _conf, K2I, TM, true);

    uint64_t nn, tn, tm;
    rocksdb::Iterator* it;

    nn = tn = tm = 0;
    it = seq[NN]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) nn++;
    it = seq[TN]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tn++;
    it = seq[TM]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tm++;
    cout << "Size SEQ: " << nn << " " << tn << " " << tm << endl;

    nn = tn = tm = 0;
    it = k2i[NN]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) nn++;
    it = k2i[TN]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tn++;
    it = k2i[TM]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tm++;
    cout << "Size K2I: " << nn << " " << tn << " " << tm << endl;

    tm = 0;
    it = i2p->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tm++;
    cout << "Size I2P: " << tm << endl;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Time merge/stats: " << time.count() << endl;
}
