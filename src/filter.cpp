/*
 * Copyright © 2015-2018 Barcelona Supercomputing Center (BSC)
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

#include "filter.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include <unordered_map>

#include "registry.hpp"
#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

filter::filter(const sm_config &conf) : stage(conf)
{
    _input_queue = sm::input_queues.at(_conf.input_format)(conf);
    _input_queue->init();

    string name = _conf.index_format;
    if (sm::index_formats.find(name) != sm::index_formats.end()) {
        cout << "Initialize: filter-" << name << endl;
        _format = sm::index_formats.at(name)(_conf);
    } else {
        cout << "Unknown filter format: " << name << endl;
        exit(1);
    }

    _executable["run"] = std::bind(&filter::run, this);
    _executable["dump"] = std::bind(&filter::dump, this);
    _executable["stats"] = std::bind(&filter::stats, this);
}

void filter::chain(const stage* prev)
{
    cout << "Chain: filter" << endl;
    _count = static_cast<const count*>(prev);
}

void filter::run()
{
    cout << "Filter: " << _conf.num_filters << " threads x "
         << _conf.num_indexes << " indexes" << endl;
    spawn("filter", std::bind(&filter::load, this, std::placeholders::_1),
          _conf.num_filters);
}

void filter::stats()
{
    _format->stats();
}

void filter::dump()
{
    _format->dump();
}

void filter::load(int fid)
{
    sm_chunk chunk;
    while (_input_queue->len > 0) {
        while (_input_queue->try_dequeue(chunk)) {
            load_chunk(fid, chunk);
            _input_queue->len--;
        }
    }
}

void filter::load_chunk(int fid, const sm_chunk &chunk)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    input_iterator *it;
    uint64_t num_reads = 0;
    sm_read read;
    sm_bulk_msg bulks[MAX_STORERS];

    it = sm::input_iterators.at(_conf.input_format)(_conf, chunk);
    while (it->next(&read)) {
        num_reads++;

        for (int i = 0; i < read.num_splits; i++) {
            int p = read.splits[i][0];
            int n = read.splits[i][1];
            if (chunk.kind == CANCER_READ)
                filter_cancer(fid, &read, &read.seq[p], n);
            else
                filter_normal(fid, &read, &read.seq[p], n);
        }

        if (num_reads % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "F: " << fid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }

        if (num_reads % 10000000 == 0) {
            bool f = _format->flush();
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "W: " << fid << " " << time.count() << " " << f << endl;
            start = std::chrono::system_clock::now();
        }
    }
}

void filter::filter_normal(int fid, const sm_read *read, const char *sub,
                           int len)
{
    if (len < _conf.k)
        return;

    int stem_len = _conf.k - 2;
    char stem[stem_len + 1];
    char kmer[_conf.k + 1];
    for (int i = 0; i <= len - _conf.k; i++) {
        strncpy(stem, &sub[i + 1], stem_len);
        stem[stem_len] = '\0';
        int order = min_order(stem, stem_len);
        if (order)
            revcomp(stem, stem_len);

        sm_table::const_iterator it;
        if (get_value(stem, &it) != 0)
            continue;

        strncpy(kmer, &sub[i], _conf.k);
        kmer[_conf.k] = '\0';
        filter_all(fid, read, i, kmer, DIR_A, order, it->second, NN);

        revcomp(kmer, _conf.k);
        order = (order + 1) % 2;
        filter_all(fid, read, i, kmer, DIR_B, order, it->second, NN);
    }
}

void filter::filter_cancer(int fid, const sm_read *read, const char *sub,
                           int len)
{
    if (len < _conf.k)
        return;

    int stem_len = _conf.k - 2;
    char stem[stem_len + 1];
    char kmer[_conf.k + 1];
    for (int i = 0; i <= len - _conf.k; i++) {
        strncpy(stem, &sub[i + 1], stem_len);
        stem[stem_len] = '\0';
        int order = min_order(stem, stem_len);
        if (order)
            revcomp(stem, stem_len);

        sm_table::const_iterator it;
        if (get_value(stem, &it) != 0)
            continue;

        strncpy(kmer, &sub[i], _conf.k);
        kmer[_conf.k] = '\0';
        filter_branch(fid, read, i, kmer, DIR_A, order, it->second, TM);
        filter_all(fid, read, i, kmer, DIR_A, order, it->second, TN);

        revcomp(kmer, _conf.k);
        order = (order + 1) % 2;
        filter_branch(fid, read, i, kmer, DIR_B, order, it->second, TM);
        filter_all(fid, read, i, kmer, DIR_B, order, it->second, TN);
    }
}

int filter::get_value(char stem[], sm_table::const_iterator *it)
{
    uint64_t m = 0;
    memcpy(&m, stem, MAP_LEN);
    map_mer(m);

    if (map_l1[m] != _conf.pid)
        return -1;
    int sid = map_l2[m];
    sm_key key = strtob4(stem);

    const sm_table* table = (*_count)[sid];
    *it = table->find(key);
    if (*it == table->end())
        return -1;
    return 0;
}

void filter::filter_all(int fid, const sm_read *read, int pos, char kmer[],
                        sm_dir dir, int order, const sm_value &counts,
                        sm_idx_set set)
{
    char first = kmer[0];
    char last = kmer[_conf.k - 1];

    for (int f = 0; f < 4; f++) {
        kmer[0] = sm::alpha[f];
        uint32_t nsum = 0;
        uint32_t tsum = 0;
        for (int l = 0; l < 4; l++) {
            nsum += counts.v[order][f][l][NORMAL_READ];
            tsum += counts.v[order][f][l][CANCER_READ];
        }
        for (int l = 0; l < 4; l++) {
            kmer[_conf.k - 1] = sm::alpha[l];
            int orderb = (order + 1) % 2;
            int fb = (3 - l);
            int lb = (3 - f);
            uint32_t na = counts.v[order ][f ][l ][NORMAL_READ];
            uint32_t ta = counts.v[order ][f ][l ][CANCER_READ];
            uint32_t nb = counts.v[orderb][fb][lb][NORMAL_READ];
            uint32_t tb = counts.v[orderb][fb][lb][CANCER_READ];
            filter_kmer(fid, read, pos, kmer, dir, na, ta, nb, tb, nsum, tsum, set);
        }
    }

    kmer[0] = first;
    kmer[_conf.k - 1] = last;
}

void filter::filter_branch(int fid, const sm_read *read, int pos, char kmer[],
                           sm_dir dir, int order, const sm_value &counts,
                           sm_idx_set set)
{
    int f = sm::code[kmer[0]] - '0';
    int l = sm::code[kmer[_conf.k - 1]] - '0';

    int orderb = (order + 1) % 2;
    int fb = (3 - l);
    int lb = (3 - f);

    uint32_t na = counts.v[order ][f ][l ][NORMAL_READ];
    uint32_t ta = counts.v[order ][f ][l ][CANCER_READ];
    uint32_t nb = counts.v[orderb][fb][lb][NORMAL_READ];
    uint32_t tb = counts.v[orderb][fb][lb][CANCER_READ];

    uint32_t nsum = 0;
    uint32_t tsum = 0;
    for (l = 0; l < 4; l++) {
        nsum += counts.v[order][f][l][NORMAL_READ];
        tsum += counts.v[order][f][l][CANCER_READ];
    }

    filter_kmer(fid, read, pos, kmer, dir, na, ta, nb, tb, nsum, tsum, set);
}

void filter::filter_kmer(int fid, const sm_read *read, int pos, char kmer[],
                         sm_dir dir, uint32_t na, uint32_t ta, uint32_t nb,
                         uint32_t tb, uint32_t nsum, uint32_t tsum,
                         sm_idx_set set)
{
    if (ta >= _conf.min_tc_a && na <= _conf.max_nc_a &&
        tb >= _conf.min_tc_b && nb <= _conf.max_nc_b) {
        if (dir == DIR_B) {
            // Recalculate reverse-complement position since the loops, and
            // thus the passed `pos', follow the forward sequence.
            pos = read->len - _conf.k - pos;
        }
        _format->update(fid, read, pos, kmer, dir, set);
    }
}
