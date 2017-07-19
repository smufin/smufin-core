#include "merge.hpp"

#include <stdio.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>

#include "index_iterator.hpp"
#include "registry.hpp"
#include "util.hpp"

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
    _executable["run_i2p_tm"] = std::bind(&merge::load, this, I2P, TM);

    _executable["stats"] = std::bind(&merge::stats, this);
    _executable["to_fastq"] = std::bind(&merge::to_fastq, this);
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
    for (int i = 0; i < _conf.num_partitions; i++) {
        for (int j = 0; j < _conf.num_indexes; j++) {
            _index_queue.enqueue(std::pair<int, int>(i, j));
        }
    }

    std::map<sm_idx_type, load_f> load_map;
    load_map[SEQ] = std::bind(&merge::load_seq, this, _1, _2, _3, _4);
    load_map[K2I] = std::bind(&merge::load_k2i, this, _1, _2, _3, _4);
    load_map[I2P] = std::bind(&merge::load_i2p, this, _1, _2, _3, _4);

    rdb_handle rdb;
    open_index_full_load(_conf, type, set, rdb);

    std::vector<std::thread> threads;
    for (int i = 0; i < _conf.num_mergers; i++)
        threads.push_back(std::thread(std::bind(&merge::load_indexes, this,
                          load_map[type], rdb, set)));
    cout << "Spawned " << threads.size() << " merger threads" << endl;
    for (auto& thread: threads)
        thread.join();

    delete rdb.cfs[0];
    delete rdb.db;
}

// Load all partitions of a certain type using `load_part' and the given set
// to a RocksDB database.
void merge::load_indexes(load_f load_index, rdb_handle &rdb, sm_idx_set set)
{
    std::pair<int, int> id;
    while (_index_queue.try_dequeue(id)) {
        load_index(rdb, set, id.first, id.second);
    }
}

// Load SEQ index data for a given set and partition `pid' to the database.
void merge::load_seq(rdb_handle &rdb, sm_idx_set set, int pid, int iid)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    index_iterator<seq_t>* it;
    it = sm::seq_iterators.at(_conf.index_format)(_conf, set, pid, iid);
    if (!it->init())
        return;

    uint64_t n = 0;
    rocksdb::WriteBatch batch;
    rocksdb::WriteOptions w_options;
    w_options.disableWAL = true;

    while (it->next()) {
        const seq_t *i = it->get();
        batch.Put(rdb.cfs[0], i->first, i->second);
        if (n % 10000 == 0) {
            rdb.db->Write(w_options, &batch);
            batch.Clear();
        }
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << pid << " " << iid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }

    rdb.db->Write(w_options, &batch);
}

// Load K2I index data for a given set and partition `pid' to the database.
void merge::load_k2i(rdb_handle &rdb, sm_idx_set set, int pid, int iid)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    index_iterator<k2i_t>* it;
    it = sm::k2i_iterators.at(_conf.index_format)(_conf, set, pid, iid);
    if (!it->init())
        return;

    uint64_t n = 0;
    rocksdb::WriteOptions w_options;
    w_options.disableWAL = true;

    while (it->next()) {
        const k2i_t *i = it->get();
        rdb.db->Merge(w_options, rdb.cfs[0], i->first, i->second);
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << pid << " " << iid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
}

// Load I2P index data for a given set and partition `pid' to the database.
void merge::load_i2p(rdb_handle &rdb, sm_idx_set set, int pid, int iid)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    index_iterator<i2p_t>* it;
    it = sm::i2p_iterators.at(_conf.index_format)(_conf, set, pid, iid);
    if (!it->init())
        return;

    uint64_t n = 0;
    rocksdb::WriteOptions w_options;
    w_options.disableWAL = true;

    while (it->next()) {
        const i2p_t *i = it->get();
        string serialized;
        encode_pos(i->second, serialized);
        rdb.db->Merge(w_options, rdb.cfs[0], i->first, serialized);
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << pid << " " << iid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
}

void merge::stats()
{
    rdb_handle seq[NUM_SETS];
    rdb_handle k2i[2];
    rdb_handle i2p;

    open_index_full_iter(_conf, SEQ, NN, seq[NN]);
    open_index_full_iter(_conf, SEQ, TN, seq[TN]);
    open_index_full_iter(_conf, SEQ, TM, seq[TM]);
    open_index_full_iter(_conf, K2I, NN, k2i[NN]);
    open_index_full_iter(_conf, K2I, TN, k2i[TN]);
    open_index_full_iter(_conf, I2P, TM, i2p);

    rocksdb::ReadOptions r_opt;
    uint64_t nn, tn, tm;
    rocksdb::Iterator* it;

    nn = tn = tm = 0;
    it = seq[NN].db->NewIterator(r_opt, seq[NN].cfs[0]);
    for (it->SeekToFirst(); it->Valid(); it->Next()) nn++;
    it = seq[TN].db->NewIterator(r_opt, seq[TN].cfs[0]);
    for (it->SeekToFirst(); it->Valid(); it->Next()) tn++;
    it = seq[TM].db->NewIterator(r_opt, seq[TM].cfs[0]);
    for (it->SeekToFirst(); it->Valid(); it->Next()) tm++;
    cout << "Size SEQ: " << nn << " " << tn << " " << tm << endl;

    nn = tn = 0;
    it = k2i[NN].db->NewIterator(r_opt, k2i[NN].cfs[0]);
    for (it->SeekToFirst(); it->Valid(); it->Next()) nn++;
    it = k2i[TN].db->NewIterator(r_opt, k2i[TN].cfs[0]);
    for (it->SeekToFirst(); it->Valid(); it->Next()) tn++;
    cout << "Size K2I: " << nn << " " << tn << endl;

    tm = 0;
    it = i2p.db->NewIterator(r_opt, i2p.cfs[0]);
    for (it->SeekToFirst(); it->Valid(); it->Next()) tm++;
    cout << "Size I2P: " << tm << endl;
}

void merge::to_fastq()
{
    std::vector<sm_idx_set> list = { NN, TN };
    spawn<sm_idx_set>("to_fastq", std::bind(&merge::to_fastq_set, this,
                      std::placeholders::_1), list);
}

void merge::to_fastq_set(sm_idx_set set)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path << "/merge-" << sm::sets[set] << ".fastq";
    ofs.open(file.str());

    rdb_handle rdb;
    open_index_full_iter(_conf, SEQ, set, rdb);

    rocksdb::ReadOptions r_opt;
    rocksdb::Iterator* it = rdb.db->NewIterator(r_opt, rdb.cfs[0]);
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        string id = it->key().ToString();
        string seq = it->value().ToString();
        // Quality is not available at this point, so score set to lowest
        // quality and should be ignored.
        string qual = string(seq.length(), '!');
        ofs << "@" << id << "\n" << seq << "\n+\n" << qual << "\n";
    }

    delete rdb.cfs[0];
    delete rdb.db;
    ofs.close();
}
