/*
 * Copyright © 2015-2019 Barcelona Supercomputing Center (BSC)
 *
 * This file is part of SMUFIN Core. SMUFIN Core is released under the SMUFIN
 * Public License, and may not be used except in compliance with it. This file
 * is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; see the SMUFIN Public License for more details. You should have
 * received a copy of the SMUFIN Public License along with this file. If not,
 * see <https://github.com/smufin/smufin-core/blob/master/COPYING>.
 *
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2018
 */

#include "prune.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "registry.hpp"
#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

prune::prune(const sm_config &conf) : stage(conf)
{
    _input_queue = sm::input_queues.at(_conf.input_format)(conf);
    _input_queue->init(_conf.num_loaders);

    _all_size = _conf.all_size / _conf.num_partitions / _conf.num_storers;
    _allowed_size = _conf.allowed_size / _conf.num_partitions /
                    _conf.num_storers;

    _executable["run"] = std::bind(&prune::run, this);
}

prune::~prune()
{
    for (int i = 0; i < _conf.num_loaders; i++) {
        delete _allowed[i];
        for (int j = 0; j < _conf.num_loaders; j++) {
            delete _queues[i][j];
        }
    }
}

void prune::run()
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
            _queues[i][j] = new sm_prune_queue(PRUNE_QUEUE_LEN);
        }
    }

    std::vector<std::thread> loaders;
    for (int i = 0; i < _conf.num_loaders; i++)
        loaders.push_back(std::thread(&prune::load, this, i));
    cout << "Spawned " << loaders.size() << " prune loader threads" << endl;

    std::vector<std::thread> storers;
    for (int i = 0; i < _conf.num_storers; i++)
        storers.push_back(std::thread(&prune::add, this, i));
    cout << "Spawned " << storers.size() << " prune storer threads" << endl;

    for (auto& loader: loaders)
        loader.join();
    _done = true;
    for (auto& storer: storers)
        storer.join();
}

void prune::load(int lid)
{
    sm_chunk chunk;
    while (_input_queue->len > 0) {
        while (_input_queue->try_dequeue(chunk)) {
            load_chunk(lid, chunk);
            _input_queue->len--;
        }
    }
}

void prune::load_chunk(int lid, const sm_chunk &chunk)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    input_iterator *it;
    uint64_t num_reads = 0;
    sm_read read;
    sm_bulk_key bulks[MAX_STORERS];

    it = sm::input_iterators.at(_conf.input_format)(_conf, chunk);
    while (it->next(&read)) {
        num_reads++;

        for (int i = 0; i < read.num_splits; i++) {
            int p = read.splits[i][0];
            int n = read.splits[i][1];
            load_sub(lid, &read.seq[p], n, bulks);
        }

        if (num_reads % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "P: " << lid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
    }

    for (int sid = 0; sid < _conf.num_storers; sid++) {
        while (!_queues[sid][lid]->try_enqueue(bulks[sid])) {
            continue;
        }
        bulks[sid].num = 0;
    }
}

inline void prune::load_sub(int lid, const char* sub, int len,
                            sm_bulk_key* bulks)
{
    if (len < _conf.k)
        return;

    char stem_str[_conf.stem_len + 1];
    for (int i = 0; i <= len - _conf.k; i++) {
        strncpy(stem_str, &sub[i + 1], _conf.stem_len);
        stem_str[_conf.stem_len] = '\0';

        uint64_t m = 0;
        memcpy(&m, &stem_str[_conf.map_pos], MAP_LEN);
        map_mer(m);

        if (map_l1[m] != _conf.pid)
            continue;
        int sid = map_l2[m];
        sm_key stem_key = strtob4(stem_str);

        bulks[sid].array[bulks[sid].num] = stem_key;
        bulks[sid].num++;

        if (bulks[sid].num == BULK_KEY_LEN) {
            while (!_queues[sid][lid]->try_enqueue(bulks[sid])) {
                continue;
            }
            bulks[sid].num = 0;
        }
    }
}

void prune::add(int sid)
{
    double fp = _conf.false_positive_rate;
    _all[sid] = new bf::basic_bloom_filter(fp, _all_size);
    _allowed[sid] = new bf::basic_bloom_filter(fp, _allowed_size);

    sm_bulk_key* pmsg;
    while (!_done) {
        for (int lid = 0; lid < _conf.num_loaders; lid++) {
            pmsg = _queues[sid][lid]->peek();
            while (pmsg) {
                for (int i = 0; i < pmsg->num; i++) {
                    add_key(sid, pmsg->array[i]);
                }
                _queues[sid][lid]->pop();
                pmsg = _queues[sid][lid]->peek();
            }
        }
    }

    for (int lid = 0; lid < _conf.num_loaders; lid++) {
        pmsg = _queues[sid][lid]->peek();
        while (pmsg) {
            for (int i = 0; i < pmsg->num; i++) {
                add_key(sid, pmsg->array[i]);
            }
            _queues[sid][lid]->pop();
            pmsg = _queues[sid][lid]->peek();
        }
    }

    delete _all[sid];
}

void prune::add_key(int sid, sm_key stem)
{
    sm_key root = to_root(stem, _conf.stem_len);
    if (!_all[sid]->lookup(root)) {
        _all[sid]->add(root);
    } else if (!_allowed[sid]->lookup(root)) {
        _allowed[sid]->add(root);
    }
}
