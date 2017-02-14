#include "count.hpp"

#include <errno.h>
#include <math.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <thread>

#include <boost/algorithm/string.hpp>

#include "input.hpp"
#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

count::count(const sm_config &conf) : stage(conf)
{
    _input_queue = new input_queue(conf);
    _input_queue->init();

    _table_size = _conf.table_size / _conf.num_partitions / _conf.num_storers;
    _cache_size = _conf.cache_size / _conf.num_partitions / _conf.num_storers;

    _executable["run"] = std::bind(&count::run, this);
    _executable["dump"] = std::bind(&count::dump, this);
    _executable["restore"] = std::bind(&count::restore, this);
    _executable["stats"] = std::bind(&count::stats, this);
}

void count::chain(const stage* prev)
{
    cout << "Chain: count" << endl;
    _prune = static_cast<const prune*>(prev);
    _enable_prune = true;
}

void count::run()
{
    if (_conf.num_loaders > MAX_LOADERS) {
        cout << "Number of loaders is larger than MAX_LOADERS" << endl;
        exit(1);
    }

    if (_conf.num_storers > MAX_STORERS) {
        cout << "Number of storers is larger than MAX_STORERS" << endl;
        exit(1);
    }

    // Initialize message queues
    for (int i = 0; i < _conf.num_storers; i++) {
        for (int j = 0; j < _conf.num_loaders; j++) {
            _queues[i][j] = new sm_queue(COUNT_QUEUE_LEN);
        }
    }

    float table_mem = estimate_sparse(_conf.num_storers * _table_size,
                                      sizeof(sm_key), sizeof(sm_value));
    float cache_mem = estimate_sparse(_conf.num_storers * _cache_size,
                                      sizeof(sm_key), sizeof(uint8_t));
    cout << "Tables: " << _table_size << " x " << _conf.num_storers
         << " (estimated up to ~" << table_mem << "GB)" << endl;
    cout << "Caches: " << _cache_size << " x " << _conf.num_storers
         << " (estimated up to ~" << cache_mem << "GB)" << endl;

    std::vector<std::thread> loaders;
    for (int i = 0; i < _conf.num_loaders; i++)
        loaders.push_back(std::thread(&count::load, this, i));
    cout << "Spawned " << loaders.size() << " loader threads" << endl;

    std::vector<std::thread> storers;
    for (int i = 0; i < _conf.num_storers; i++)
        storers.push_back(std::thread(&count::incr, this, i));
    cout << "Spawned " << storers.size() << " storer threads" << endl;

    for (auto& loader: loaders)
        loader.join();
    _done = true;
    for (auto& storer: storers)
        storer.join();

    for (int i = 0; i < _conf.num_storers; i++)
        cout << "Table " << i << ": " << _tables[i]->size() << endl;
}

void count::load(int lid)
{
    sm_chunk chunk;
    while (_input_queue->len > 0) {
        while (_input_queue->try_dequeue(chunk)) {
            load_chunk(lid, chunk);
            _input_queue->len--;
        }
    }
}

void count::load_chunk(int lid, const sm_chunk &chunk)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    uint64_t num_reads = 0;
    sm_read read;
    sm_bulk_msg bulks[MAX_STORERS];

    input_iterator_fastq it(_conf, chunk);
    while (it.next(&read)) {
        num_reads++;

        for (int i = 0; i < read.num_splits; i++) {
            int p = read.splits[i][0];
            int n = read.splits[i][1];
            load_sub(lid, &read.seq[p], n, chunk.kind, bulks);
        }

        if (num_reads % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "C: " << lid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
    }

    for (int sid = 0; sid < _conf.num_storers; sid++) {
        while (!_queues[sid][lid]->write(bulks[sid])) {
            continue;
        }
        bulks[sid].num = 0;
    }
}

inline void count::load_sub(int lid, const char* sub, int len,
                            sm_read_kind kind, sm_bulk_msg* bulks)
{
    if (len < _conf.k)
        return;

    int stem_len = _conf.k - 2;
    char stem[stem_len + 1];
    for (int i = 0; i <= len - _conf.k; i++) {
        strncpy(stem, &sub[i + 1], stem_len);
        stem[stem_len] = '\0';

        uint64_t m = 0;
        memcpy(&m, stem, MAP_LEN);
        hash_5mer(m);

        if (map_l1[m] != _conf.pid)
            continue;
        int sid = map_l2[m];
        sm_key key = strtob4(stem);

        if (_enable_prune && !(*_prune)[sid]->lookup(key)) {
            continue;
        }

        sm_value_offset off;
        off.first = sm::code[sub[i]] - '0';
        off.last = sm::code[sub[i + _conf.k - 1]] - '0';
        off.kind = kind;

        bulks[sid].array[bulks[sid].num] = sm_msg(key, off);
        bulks[sid].num++;

        if (bulks[sid].num == BULK_MSG_LEN) {
            while (!_queues[sid][lid]->write(bulks[sid])) {
                continue;
            }
            bulks[sid].num = 0;
        }
    }
}

void count::incr(int sid)
{
    _tables[sid] = new sm_table();
    _tables[sid]->resize(_table_size);

    if (_conf.enable_cache) {
        _caches[sid] = new sm_cache();
        _caches[sid]->resize(_cache_size);
    }

    sm_bulk_msg* pmsg;
    while (!_done) {
        for (int lid = 0; lid < _conf.num_loaders; lid++) {
            pmsg = _queues[sid][lid]->frontPtr();
            while (pmsg) {
                for (int i = 0; i < pmsg->num; i++) {
                    incr_key(sid, pmsg->array[i].first, pmsg->array[i].second);
                }
                _queues[sid][lid]->popFront();
                pmsg = _queues[sid][lid]->frontPtr();
            }
        }
    }

    for (int lid = 0; lid < _conf.num_loaders; lid++) {
        pmsg = _queues[sid][lid]->frontPtr();
        while (pmsg) {
            for (int i = 0; i < pmsg->num; i++) {
                incr_key(sid, pmsg->array[i].first, pmsg->array[i].second);
            }
            _queues[sid][lid]->popFront();
            pmsg = _queues[sid][lid]->frontPtr();
        }
    }

    if (_conf.enable_cache) {
        delete _caches[sid];
    }
}

inline void count::incr_key(int sid, sm_key key, sm_value_offset off)
{
    // Use sm_cache to hold keys with a single appearance; as soon as a key in
    // increased more than once, it is placed into sm_table. The steps are as
    // follows:
    //
    // - Find key in cache.
    //   - Key doesn't exist in cache: insert in cache.
    //   - Key exists in cache: find key in table.
    //     - Key doesn't exist in table: insert key and cache value in table.
    //     - Key exists in table: update entry if there's no overflow.

    sm_cache::const_iterator cit;

    if (_conf.enable_cache) {
        cit = _caches[sid]->find(key);
        if (cit == _caches[sid]->end()) {
            uint8_t val = (off.first << 6) | (off.last << 4) | (off.kind << 2);
            _caches[sid]->insert(std::pair<sm_key, uint8_t>(key, val));
            return;
        }
    }

    sm_table::const_iterator it = _tables[sid]->find(key);
    if (it == _tables[sid]->end()) {
        sm_value val;
        if (_conf.enable_cache) {
            uint8_t cache_value = cit->second;
            sm_value_offset coff;
            coff.first = cache_value >> 6;
            coff.last = (cache_value >> 4) & 0x03;
            coff.kind = (sm_read_kind) ((cache_value >> 2) & 0x03);
            val.v[coff.first][coff.last][coff.kind] = 1;
        }
        val.v[off.first][off.last][off.kind]++;
        _tables[sid]->insert(std::pair<sm_key, sm_value>(key, val));
    } else {
        uint32_t inc = it->second.v[off.first][off.last][off.kind] + 1;
        uint16_t over = inc >> 16;
        uint16_t count = inc & 0x0000FFFF;
        if (over == 0)
            (*_tables[sid])[key].v[off.first][off.last][off.kind] = count;
    }
}

void count::dump()
{
    spawn("dump", std::bind(&count::dump_table, this, std::placeholders::_1),
          _conf.num_storers);
}

void count::dump_table(int sid)
{
    std::ostringstream fs;
    fs << _conf.output_path << "/table." << _conf.pid << "-" << sid << ".sht";
    string file = fs.str();
    cout << "Serialize " << file << endl;

    FILE* fp = fopen(file.c_str(), "w");
    if (fp == NULL) {
        cout << "Failed to open " << file << " (" << errno << ")" << endl;
        exit(1);
    }

    if (!_tables[sid]->serialize(sm_table::NopointerSerializer(), fp)) {
        cout << "Failed to serialize table " << _conf.pid << "-" << sid << endl;
        exit(1);
    }

    fclose(fp);
}

void count::restore()
{
    spawn("restore", std::bind(&count::restore_table, this, std::placeholders::_1),
          _conf.num_storers);
}

void count::restore_table(int sid)
{
    _tables[sid] = new sm_table();
    _tables[sid]->resize(_table_size);

    std::ostringstream fs;
    fs << _conf.output_path << "/table." << _conf.pid << "-" << sid << ".sht";
    string file = fs.str();
    cout << "Unserialize " << file << endl;

    FILE* fp = fopen(file.c_str(), "r");
    if (fp == NULL) {
        cout << "Failed to open " << file << " (" << errno << ")" << endl;
        exit(0);
    }

    _tables[sid]->unserialize(sm_table::NopointerSerializer(), fp);
    fclose(fp);
}

void count::stats()
{
    std::map<uint64_t, uint64_t> hist_n, hist_t;

    uint64_t total_stems = 0;
    uint64_t total_once = 0;
    uint64_t total_kmers = 0;
    uint64_t total_sum = 0;

    for (int i = 0; i < _conf.num_storers; i++) {
        uint64_t num_stems = _tables[i]->size();
        uint64_t num_once = 0;
        uint64_t num_kmers = 0;
        uint64_t sum = 0;
        for (sm_table::const_iterator it = _tables[i]->begin();
             it != _tables[i]->end(); ++it) {
            uint64_t sum_t = 0;
            uint64_t sum_n = 0;
            for (int f = 0; f < 4; f++) {
                for (int l = 0; l < 4; l++) {
                    uint16_t nc = it->second.v[f][l][NORMAL_READ];
                    uint16_t tc = it->second.v[f][l][CANCER_READ];
                    sum_n += nc;
                    sum_t += tc;
                    num_kmers += (nc + tc > 0) ? 1 : 0;
                }
            }
            if (sum_n > 0)
                hist_n[int(log2(sum_n))]++;
            if (sum_t > 0)
                hist_t[int(log2(sum_t))]++;
            if (sum_n + sum_t == 1)
                num_once++;
            sum += sum_n + sum_t;
        }
        cout << "Table " << i << ": " << num_stems << " " << num_once << " "
             << num_kmers << " " << sum << endl;
        total_stems += num_stems;
        total_once += num_once;
        total_kmers += num_kmers;
        total_sum += sum;
    }

    for (auto& h: hist_n) {
        cout << "Histo N: " << h.first << " " << exp2(h.first)
             << " " << h.second << endl;
    }

    for (auto& h: hist_t) {
        cout << "Histo T: " << h.first << " " << exp2(h.first)
             << " " << h.second << endl;
    }

    cout << "Number of stems: " << total_stems << endl;
    cout << "Number of stems seen once: " << total_once << endl;
    cout << "Number of kmers: " << total_kmers << endl;
    cout << "Sum of counters: " << total_sum << endl;
}
