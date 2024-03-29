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
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2019
 */

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

#include "filter.hpp"
#include "registry.hpp"
#include "input_iterator_fastq.hpp"
#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

count::count(const sm_config &conf) : stage(conf)
{
    _input_queue = sm::input_queues.at(_conf.input_format)(conf);
    _input_queue->init(_conf.num_loaders);

    _table_size = _conf.table_size / _conf.num_partitions / _conf.num_storers;
    _cache_size = _conf.cache_size / _conf.num_partitions / _conf.num_storers;

    _executable["run"] = std::bind(&count::run, this);
    _executable["dump"] = std::bind(&count::dump, this);
    _executable["restore"] = std::bind(&count::restore, this);
    _executable["stats"] = std::bind(&count::stats, this);
    _executable["export"] = std::bind(&count::export_csv, this);
    _executable["annotate"] = std::bind(&count::annotate, this);
}

void count::chain(const stage* prev)
{
    cout << "Chain: count" << endl;
    _prune = static_cast<const prune*>(prev);
    _enable_prune = true;
}

void count::run()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;

    start = std::chrono::system_clock::now();

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
                                      sizeof(sm_key), sizeof(sm_stem));
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

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Time count/run/load: " << time.count() << endl;

    if (_enable_prune) {
        delete _prune;
    }

    start = std::chrono::system_clock::now();

    std::vector<std::thread> converters;
    for (int i = 0; i < _conf.max_conversions; i++)
        converters.push_back(std::thread(&count::convert, this));
    cout << "Spawned " << converters.size() << " convert threads ("
         << _conf.conversion_mode << ")" << endl;
    for (auto& converter: converters)
        converter.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Time count/run/convert: " << time.count() << endl;
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
        while (!_queues[sid][lid]->try_enqueue(bulks[sid])) {
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

        sm_stem_offset off;
        off.first = sm::code[sub[i]] - '0';
        off.last = sm::code[sub[i + _conf.k - 1]] - '0';
        off.kind = kind;

        bulks[sid].array[bulks[sid].num] = sm_msg(stem_key, off);
        bulks[sid].num++;

        if (bulks[sid].num == BULK_MSG_LEN) {
            while (!_queues[sid][lid]->try_enqueue(bulks[sid])) {
                continue;
            }
            bulks[sid].num = 0;
        }
    }
}

void count::incr(int sid)
{
    _stem_tables[sid] = new sm_stem_table();
    _stem_tables[sid]->resize(_table_size);

    if (_conf.slice) {
        _stem_tables[sid]->set_resizing_parameters(0.01, 0.99);
    }

    _slices[sid] = new std::vector<int>();

    if (_conf.enable_cache) {
        _root_caches[sid] = new sm_cache();
        _root_caches[sid]->resize(_cache_size);
    }

    sm_bulk_msg *pmsg;
    while (!_done) {
        for (int lid = 0; lid < _conf.num_loaders; lid++) {
            pmsg = _queues[sid][lid]->peek();
            while (pmsg) {
                for (int i = 0; i < pmsg->num; i++) {
                    incr_key(sid, pmsg->array[i].first, pmsg->array[i].second);
                }
                _queues[sid][lid]->pop();
                pmsg = _queues[sid][lid]->peek();
            }
        }

        if (_conf.slice && _stem_tables[sid]->size() > _table_size * 0.8) {
            dump_slice(sid);
            delete _stem_tables[sid];
            _stem_tables[sid] = new sm_stem_table();
            _stem_tables[sid]->resize(_table_size);
            _stem_tables[sid]->set_resizing_parameters(0.01, 0.99);
        }
    }

    for (int lid = 0; lid < _conf.num_loaders; lid++) {
        pmsg = _queues[sid][lid]->peek();
        while (pmsg) {
            for (int i = 0; i < pmsg->num; i++) {
                incr_key(sid, pmsg->array[i].first, pmsg->array[i].second);
            }
            _queues[sid][lid]->pop();
            pmsg = _queues[sid][lid]->peek();
        }
    }

    if (_conf.enable_cache) {
        delete _root_caches[sid];
    }
}

inline void count::incr_key(int sid, sm_key stem, sm_stem_offset off)
{
    // Use sm_cache to hold keys with a single appearance; as soon as a key in
    // increased more than once, it is placed into sm_table. The steps are as
    // follows:
    //
    // - Find the stem's root in cache.
    //   - Stem's root doesn't exist in cache: insert in cache.
    //   - Stem's root exists in cache: find stem in table.
    //     - Stem doesn't exist in table:
    //       - Insert cache value in table.
    //       - Insert stem in table.
    //     - Stem exists in table: update entry if there's no overflow.

    sm_cache::const_iterator cit;

    int order = 0;
    sm_key root = to_root(stem, _conf.stem_len);
    if (root < stem) {
        order = 1;
    }

    if (_enable_prune && !(*_prune)[sid]->lookup(root)) {
        return;
    }

    if (_conf.enable_cache) {
        cit = _root_caches[sid]->find(root);
        if (cit == _root_caches[sid]->end()) {
            uint8_t val = (order << 6) | (off.first << 4) | (off.last << 2) |
                          off.kind;
            _root_caches[sid]->insert(std::pair<sm_key, uint8_t>(root, val));
            return;
        }
    }

    sm_stem_table::const_iterator it = _stem_tables[sid]->find(stem);
    if (it == _stem_tables[sid]->end()) {
        sm_stem val;

        if (_conf.enable_cache) {
            int corder = 0;
            uint8_t cache_value = cit->second;
            sm_stem_offset coff;
            uint8_t saved = (cache_value >> 7) & 0x01;
            corder = (cache_value >> 6) & 0x01;
            coff.first = (cache_value >> 4) & 0x03;
            coff.last = (cache_value >> 2) & 0x03;
            coff.kind = (sm_read_kind) (cache_value & 0x03);

            // Insert in table the kmer stored in the cache. When the stem in
            // cache is different than that of the current stem, an additional
            // insertion is needed; otherwise just insert along with current
            // stem.
            if (saved == 0) {
                if (order != corder) {
                    sm_key cstem = revcomp_code(stem, _conf.stem_len);
                    sm_stem cval;
                    cval.v[coff.first][coff.last][coff.kind] = 1;
                    _stem_tables[sid]->insert(std::pair<sm_key, sm_stem>(cstem, cval));
                } else {
                    val.v[coff.first][coff.last][coff.kind] = 1;
                }

                cache_value |= (1 << 7);
                (*_root_caches[sid])[root] = cache_value;
            }
        }

        val.v[off.first][off.last][off.kind]++;
        _stem_tables[sid]->insert(std::pair<sm_key, sm_stem>(stem, val));
    } else {
        uint32_t inc = it->second.v[off.first][off.last][off.kind] + 1;
        uint16_t over = inc >> 16;
        uint16_t count = inc & 0x0000FFFF;
        if (over == 0)
            (*_stem_tables[sid])[stem].v[off.first][off.last][off.kind] = count;
    }
}

void count::convert()
{
    while (_convert < _conf.num_storers) {
        int sid = _convert;
        bool inc = _convert.compare_exchange_weak(sid, sid + 1);
        if (sid < _conf.num_storers && inc) {
            if (_conf.conversion_mode == "mem") {
                convert_table_mem(sid);
            } else {
                convert_table_slice(sid);
            }
        }
    }
}

// In-memory conversion of a stem-indexed table to a root-indexed one.
void count::convert_table_mem(int sid)
{
    _root_tables[sid] = new sm_root_table();

    uint64_t num_stems_pass = 0;

    for (const auto& stem: *_stem_tables[sid]) {
        int order = 0;
        sm_key root = to_root(stem.first, _conf.stem_len);
        if (root < stem.first) {
            order = 1;
        }

        if (!_conf.prefilter || filter::filter_stem(_conf, stem.second)) {
            (*_root_tables[sid])[root].s[order] = stem.second;
            num_stems_pass++;
        }
    }

    uint64_t num_stems = _stem_tables[sid]->size();
    uint64_t num_roots = _root_tables[sid]->size();

    delete _stem_tables[sid];

    if (_conf.prefilter) {
        prefilter_table(sid);
    }

    cout << "Convert " << sid << ": " << num_stems << " " << num_stems_pass
         << " " << num_roots << " " << _root_tables[sid]->size() << endl;
}

// Convert stem table serializing it to disk and then reading sequentially
// without loading the whole stem table to memory.
void count::convert_table_slice(int sid)
{
    dump_slice(sid);
    delete _stem_tables[sid];

    _root_tables[sid] = new sm_root_table();

    uint64_t num_stems = 0;
    uint64_t num_stems_pass = 0;
    const bool prefilter_stem = _conf.prefilter && !_conf.slice;

    for (int i = 0; i < _slices[sid]->size(); i++) {
        std::ostringstream fs;
        fs << _conf.output_path_count << "/slice." << _conf.pid << "-" << sid
           << "." << i << ".sht";
        string file = fs.str();
        cout << "Reading " << file << " (" << (*_slices[sid])[i] << " stems)"
             << endl;

        FILE* fp = fopen(file.c_str(), "r");
        if (fp == NULL) {
            cout << "Failed to open " << file << " (" << errno << ")" << endl;
            exit(0);
        }

        uint64_t magic_num = 0;
        uint64_t table_size = 0;
        uint64_t num_buckets = 0;

        read_be(fp, &magic_num);
        read_be(fp, &table_size);
        read_be(fp, &num_buckets);

        if (magic_num != 0x24687531) {
            cout << "Failed to read " << file << " (version mismatch)" << endl;
            exit(0);
        }

        // Skip metadata: each group contains the number of non-empty buckets
        // (first 16 bits), and a bitmap (remaining 48 bits). This information is
        // ignored for now since we are mainly interested in the content of the
        // table.
        uint64_t num_groups = (table_size == 0) ? 0 : ((table_size - 1) / 48) + 1;
        fseek(fp, num_groups * 8, SEEK_CUR);

        // Read data: each stem key-value pair is loaded into the root table.
        std::pair<sm_key, sm_stem> stem;
        for (uint64_t i = 0; i < table_size; i++) {
            if (!fread(&stem, sizeof(stem), 1, fp)) {
                break;
            }

            num_stems++;

            int order = 0;
            sm_key root = to_root(stem.first, _conf.stem_len);
            if (root < stem.first) {
                order = 1;
            }

            sm_root sum;
            sm_root_table::const_iterator it = _root_tables[sid]->find(root);
            if (it != _root_tables[sid]->end()) {
                sum = it->second;
                sum.s[order] = sum.s[order] + stem.second;
            } else {
                sum.s[order] = stem.second;
            }

            if (!prefilter_stem || filter::filter_stem(_conf, stem.second)) {
                (*_root_tables[sid])[root] = sum;
                num_stems_pass++;
            }
        }

        fclose(fp);
    }

    uint64_t num_roots = _root_tables[sid]->size();

    if (_conf.prefilter) {
        prefilter_table(sid);
    }

    cout << "Convert " << sid << ": " << num_stems << " " << num_stems_pass
         << " " << num_roots << " " << _root_tables[sid]->size() << endl;
}

void count::prefilter_table(int sid)
{
    sm_root_table* table = new sm_root_table();

    for (const auto& root: *_root_tables[sid]) {
        if (filter::filter_root(_conf, root.second)) {
            table->insert(root);
        }
    }

    delete _root_tables[sid];
    _root_tables[sid] = table;
}

void count::dump()
{
    spawn("dump", std::bind(&count::dump_table, this, std::placeholders::_1),
          _conf.num_storers);
}

void count::dump_table(int sid)
{
    std::ostringstream fs;
    fs << _conf.output_path_count << "/table." << _conf.pid << "-" << sid
       << ".sht";
    string file = fs.str();
    cout << "Serialize " << file << endl;

    FILE* fp = fopen(file.c_str(), "w");
    if (fp == NULL) {
        cout << "Failed to open " << file << " (" << errno << ")" << endl;
        exit(1);
    }

    if (!_root_tables[sid]->serialize(sm_root_table::NopointerSerializer(), fp)) {
        cout << "Failed to serialize table " << _conf.pid << "-" << sid << endl;
        exit(1);
    }

    fclose(fp);
}

void count::dump_slice(int sid)
{
    if (_stem_tables[sid]->size() == 0) {
        return;
    }

    std::ostringstream fs;
    fs << _conf.output_path_count << "/slice." << _conf.pid << "-" << sid
       << "." << _slices[sid]->size() << ".sht";
    string file = fs.str();
    cout << "Serialize " << file << endl;

    FILE* fp = fopen(file.c_str(), "w");
    if (fp == NULL) {
        cout << "Failed to open " << file << " (" << errno << ")" << endl;
        exit(1);
    }

    if (!_stem_tables[sid]->serialize(sm_stem_table::NopointerSerializer(), fp)) {
        cout << "Failed to serialize slice " << _conf.pid << "-" << sid << endl;
        exit(1);
    }

    fclose(fp);

    _slices[sid]->push_back(_stem_tables[sid]->size());
}

void count::restore()
{
    spawn("restore", std::bind(&count::restore_table, this,
          std::placeholders::_1), _conf.num_storers);
}

void count::restore_table(int sid)
{
    _root_tables[sid] = new sm_root_table();

    std::ostringstream fs;
    fs << _conf.output_path_count << "/table." << _conf.pid << "-" << sid
       << ".sht";
    string file = fs.str();
    cout << "Unserialize " << file << endl;

    FILE* fp = fopen(file.c_str(), "r");
    if (fp == NULL) {
        cout << "Failed to open " << file << " (" << errno << ")" << endl;
        exit(0);
    }

    _root_tables[sid]->unserialize(sm_root_table::NopointerSerializer(), fp);
    fclose(fp);
}

void count::stats()
{
    std::map<uint64_t, uint64_t> hist_n, hist_t;

    uint64_t total_roots = 0;
    uint64_t total_stems = 0;
    uint64_t total_once = 0;
    uint64_t total_kmers = 0;
    uint64_t total_sum = 0;

    uint64_t total_hits_roots = 0;
    uint64_t total_hits_stems = 0;
    uint64_t total_hits_kmers = 0;

    for (int i = 0; i < _conf.num_storers; i++) {
        uint64_t num_roots = _root_tables[i]->size();
        uint64_t num_stems = 0;
        uint64_t num_once = 0;
        uint64_t num_kmers = 0;
        uint64_t sum = 0;
        uint64_t num_hits_roots = 0;
        uint64_t num_hits_stems = 0;
        uint64_t num_hits_kmers = 0;
        for (sm_root_table::const_iterator it = _root_tables[i]->begin();
             it != _root_tables[i]->end(); ++it) {
            uint64_t hits[2] = {0};
            for (int order = 0; order < 2; order++) {
                uint64_t sum_t = 0;
                uint64_t sum_n = 0;
                for (int f = 0; f < 4; f++) {
                    for (int l = 0; l < 4; l++) {
                        int orderb = (order + 1) % 2;
                        int fb = sm::comp_code[l];
                        int lb = sm::comp_code[f];
                        uint16_t na = it->second.s[order ].v[f ][l ][NORMAL_READ];
                        uint16_t ta = it->second.s[order ].v[f ][l ][CANCER_READ];
                        uint16_t nb = it->second.s[orderb].v[fb][lb][NORMAL_READ];
                        uint16_t tb = it->second.s[orderb].v[fb][lb][CANCER_READ];
                        sum_n += na;
                        sum_t += ta;
                        num_kmers += (na + ta > 0) ? 1 : 0;
                        if (filter::condition(_conf, na, ta, nb, tb)) {
                            hits[order]++;
                        }
                    }
                }
                if (sum_n > 0)
                    hist_n[int(log2(sum_n))]++;
                if (sum_t > 0)
                    hist_t[int(log2(sum_t))]++;
                if (sum_n + sum_t >= 1)
                    num_stems++;
                if (sum_n + sum_t == 1)
                    num_once++;
                sum += sum_n + sum_t;

                if (hits[0] > 0 || hits[1] > 0)
                    num_hits_roots++;
                if (hits[0] > 0)
                    num_hits_stems++;
                if (hits[1] > 0)
                    num_hits_stems++;
                num_hits_kmers += hits[0] + hits[1];
            }
        }
        cout << "Table " << i << ": " << num_roots << " " << num_stems << " "
             << num_once << " " << num_kmers << " " << sum << " "
             << num_hits_roots << " " << num_hits_stems << " "
             << num_hits_kmers << endl;
        total_roots += num_roots;
        total_stems += num_stems;
        total_once += num_once;
        total_kmers += num_kmers;
        total_sum += sum;
        total_hits_roots += num_hits_roots;
        total_hits_stems += num_hits_stems;
        total_hits_kmers += num_hits_kmers;
    }

    for (auto& h: hist_n) {
        cout << "Histo N: " << h.first << " " << exp2(h.first)
             << " " << h.second << endl;
    }

    for (auto& h: hist_t) {
        cout << "Histo T: " << h.first << " " << exp2(h.first)
             << " " << h.second << endl;
    }

    cout << "Number of roots: " << total_roots << endl;
    cout << "Number of stems: " << total_stems << endl;
    cout << "Number of stems seen once: " << total_once << endl;
    cout << "Number of kmers: " << total_kmers << endl;
    cout << "Sum of counters: " << total_sum << endl;
    cout << "Number of filter hits (roots): " << total_hits_roots << endl;
    cout << "Number of filter hits (stems): " << total_hits_stems << endl;
    cout << "Number of filter hits (kmers): " << total_hits_kmers << endl;
}

void count::export_csv()
{
    spawn("export", std::bind(&count::export_csv_table, this,
          std::placeholders::_1), _conf.num_storers);
}

void count::export_csv_table(int sid)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path_count << "/table." << _conf.pid << "-" << sid
         << ".csv";
    ofs.open(file.str());

    char kmer[_conf.k + 1];
    for (const auto& root: *_root_tables[sid]) {
        for (int o = 0; o < 2; o++) {
            for (int f = 0; f < 4; f++) {
                for (int l = 0; l < 4; l++) {
                    uint16_t nc = root.second.s[o].v[f][l][NORMAL_READ];
                    uint16_t tc = root.second.s[o].v[f][l][CANCER_READ];
                    if ((nc > _conf.export_min && nc < _conf.export_max) ||
                        (tc > _conf.export_min && tc < _conf.export_max)) {
                        // Rebuild kmer based on coded stem and first/last base
                        kmer[0] = sm::alpha[f];
                        b4tostr(root.first, _conf.stem_len, &kmer[1]);
                        kmer[_conf.k - 1] = sm::alpha[l];
                        kmer[_conf.k] = '\0';
                        if (o) {
                            revcomp(&kmer[1], _conf.stem_len);
                        }
                        ofs << kmer << "," << nc << "," << tc << "\n";
                    }
                }
            }
        }
    }

    ofs.close();
}

void count::annotate()
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path_count << "/annotate." << _conf.pid << ".json";
    ofs.open(file.str());

    sm_chunk chunk;
    chunk.file = _conf.annotate_input;
    chunk.begin = -1;
    chunk.end = -1;
    chunk.kind = NORMAL_READ; // Kind not relevant in this context.

    input_iterator_fastq *it;
    sm_read read;
    it = new input_iterator_fastq(_conf, chunk);
    it->check = false;

    bool first = true;
    ofs << "{";
    while (it->next(&read)) {
        if (!first)
            ofs << ",";
        first = false;
        ofs << "\"" << read.id << "\":{";
        ofs << "\"seq\":\"" << read.seq << "\",";
        ofs << "\"kmers\":{";
        for (int i = 0; i < read.num_splits; i++) {
            int p = read.splits[i][0];
            int n = read.splits[i][1];
            annotate_sub(&read.seq[p], p, n, ofs);
        }
        ofs << "}}";
    }
    ofs << "}";
}

void count::annotate_sub(const char* sub, int pos, int len, std::ofstream &ofs)
{
    if (len < _conf.k)
        return;

    bool first = true;
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

        int order = 0;
        sm_key stem = strtob4(stem_str);
        sm_key root = to_root(stem, _conf.stem_len);
        if (root < stem) {
            order = 1;
        }

        sm_root_table::const_iterator it = _root_tables[sid]->find(root);
        if (it != _root_tables[sid]->end()) {
            uint8_t f = sm::code[sub[i]] - '0';
            uint8_t l = sm::code[sub[i + _conf.k - 1]] - '0';
            uint32_t nc = it->second.s[order].v[f][l][NORMAL_READ];
            uint32_t tc = it->second.s[order].v[f][l][CANCER_READ];

            if (!first)
                ofs << ",";
            first = false;
            ofs << "\"" << (pos + i) << "\":[" << nc << "," << tc << "]";
        }
    }
}
