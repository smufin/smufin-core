#include "filter.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include <unordered_map>

#include "filter_format_plain.hpp"
#include "filter_format_rocks.hpp"
#include "input_iterator_bam.hpp"
#include "input_iterator_fastq.hpp"
#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

filter::filter(const sm_config &conf) : stage(conf)
{
    _input_queue = new input_queue(conf);
    _input_queue->init();

    std::unordered_map<string, filter_format*(*)(const sm_config &)> formats;
    formats["plain"] = &filter_format::create<filter_format_plain>;
    formats["rocks"] = &filter_format::create<filter_format_rocks>;

    string name = _conf.filter_format;
    if (formats.find(name) != formats.end()) {
        cout << "Initialize: filter-" << name << endl;
        _format = formats[name](_conf);
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

    uint64_t num_reads = 0;
    sm_read read;
    sm_bulk_msg bulks[MAX_STORERS];

    input_iterator_fastq it(_conf, chunk);
    while (it.next(&read)) {
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

    char kmer[_conf.k + 1];
    for (int i = 0; i <= len - _conf.k; i++) {
        strncpy(kmer, &sub[i], _conf.k);
        kmer[_conf.k] = '\0';
        filter_all(fid, read, i, false, kmer, NN);

        strncpy(kmer, &sub[i], _conf.k);
        kmer[_conf.k] = '\0';
        revcomp(kmer, _conf.k);
        filter_all(fid, read, i, true, kmer, NN);
    }
}

void filter::filter_cancer(int fid, const sm_read *read, const char *sub,
                           int len)
{
    if (len < _conf.k)
        return;

    char kmer[_conf.k + 1];
    for (int i = 0; i <= len - _conf.k; i++) {
        strncpy(kmer, &sub[i], _conf.k);
        kmer[_conf.k] = '\0';
        filter_branch(fid, read, i, false, kmer, TM);
        filter_all(fid, read, i, false, kmer, TN);

        strncpy(kmer, &sub[i], _conf.k);
        kmer[_conf.k] = '\0';
        revcomp(kmer, _conf.k);
        filter_branch(fid, read, i, true, kmer, TM);
        filter_all(fid, read, i, true, kmer, TN);
    }
}

int filter::get_value(int fid, char kmer[], sm_table::const_iterator *it)
{
    char last = kmer[_conf.k - 1];
    kmer[_conf.k - 1] = '\0';

    uint64_t m = 0;
    memcpy(&m, &kmer[1], MAP_LEN);
    hash_5mer(m);

    if (map_l1[m] != _conf.pid)
        return -1;
    int sid = map_l2[m];
    sm_key key = strtob4(&kmer[1]);

    const sm_table* table = (*_count)[sid];
    *it = table->find(key);
    if (*it == table->end())
        return -1;

    kmer[_conf.k - 1] = last;
    return 0;
}

void filter::filter_branch(int fid, const sm_read *read, int pos, bool rev,
                           char kmer[], sm_idx_set set)
{
    sm_table::const_iterator it;
    if (get_value(fid, kmer, &it) != 0)
        return;
    int f = sm::code[kmer[0]] - '0';
    int l = sm::code[kmer[_conf.k - 1]] - '0';
    uint32_t nc = it->second.v[f][l][NORMAL_READ];
    uint32_t tc = it->second.v[f][l][CANCER_READ];
    uint32_t nsum = 0;
    uint32_t tsum = 0;
    for (l = 0; l < 4; l++) {
        nsum += it->second.v[f][l][NORMAL_READ];
        tsum += it->second.v[f][l][CANCER_READ];
    }
    filter_kmer(read, pos, rev, kmer, nc, tc, nsum, tsum, set);
}

void filter::filter_all(int fid, const sm_read *read, int pos, bool rev,
                        char kmer[], sm_idx_set set)
{
    sm_table::const_iterator it;
    if (get_value(fid, kmer, &it) != 0)
        return;
    for (int f = 0; f < 4; f++) {
        kmer[0] = sm::alpha[f];
        uint32_t nsum = 0;
        uint32_t tsum = 0;
        for (int l = 0; l < 4; l++) {
            nsum += it->second.v[f][l][NORMAL_READ];
            tsum += it->second.v[f][l][CANCER_READ];
        }
        for (int l = 0; l < 4; l++) {
            kmer[_conf.k - 1] = sm::alpha[l];
            uint32_t nc = it->second.v[f][l][NORMAL_READ];
            uint32_t tc = it->second.v[f][l][CANCER_READ];
            filter_kmer(read, pos, rev, kmer, nc, tc, nsum, tsum, set);
        }
    }
}

void filter::filter_kmer(const sm_read *read, int pos, bool rev, char kmer[],
                         uint32_t nc, uint32_t tc, uint32_t nsum,
                         uint32_t tsum, sm_idx_set set)
{
    if (tc >= _conf.min_tc && nc <= _conf.max_nc) {
        if (rev) {
            // Recalculate reverse-complement position since the loops, and
            // thus the passed `pos', follow the forward sequence.
            pos = read->len - _conf.k - pos;
        }
        _format->update(read, pos, rev, kmer, set);
    }
}
