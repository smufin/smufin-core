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

#ifndef __SM_COUNT_H__
#define __SM_COUNT_H__

#include <string>

#include <readerwriterqueue.h>
#include <google/sparse_hash_map>

#include "common.hpp"
#include "input.hpp"
#include "prune.hpp"
#include "stage.hpp"

#define BULK_MSG_LEN 128
#define COUNT_QUEUE_LEN 512

// A value of the hashtable that contains normal and tumoral counters for all
// inflections of a stem and its reverse complement. The multidimensional
// array `v' is indexed as follows: v[A][B][C][D], where:
//  - A stands for the direction, sorted by the encoded sequence (the stem or
//    its reverse complement) in increasing order
//  - B and C are the numeric code of the first and last bases respectively
//    (as defined by sm::code)
//  - D is the kind (as defined by sm_read_kind)
typedef struct sm_stem { uint16_t v[4][4][2] = {{{0}}}; } sm_stem;
typedef struct sm_root { sm_stem s[2]; } sm_root;

typedef google::sparse_hash_map<sm_key, uint8_t, sm_hasher<sm_key>> sm_cache;
typedef google::sparse_hash_map<sm_key, sm_stem, sm_hasher<sm_key>> sm_stem_table;
typedef google::sparse_hash_map<sm_key, sm_root, sm_hasher<sm_key>> sm_root_table;

// Contains data to calculate a position within a sm_value multidimensional
// array. `first' and `last' are integers in the range 0..3 and contain codes
// representing the first and last bases of a kmer using the same mapping as
// sm::code.
typedef struct {
    uint8_t first;
    uint8_t last;
    sm_read_kind kind;
} sm_stem_offset;

typedef std::pair<sm_key, sm_stem_offset> sm_msg;

typedef struct {
    uint16_t num = 0;
    sm_msg array[BULK_MSG_LEN];
} sm_bulk_msg;

typedef moodycamel::ReaderWriterQueue<sm_bulk_msg> sm_queue;

// Stage that reads input chunks, splits sequences into kmers, and builds a
// table of normal and tumoral kmer frequencies. `count' provides an in-memory
// implementation, and uses a cache that holds kmers that are seen only once.
// Frequency tables are indexed by stem instead of kmers, and each entry
// contains normal and tumoral counters for all inflections.
class count : public stage
{
public:
    count(const sm_config &conf);
    void run();
    void chain(const stage* prev);

    inline const sm_root_table* operator[](int sid) const {
        return _root_tables[sid];
    };

private:
    uint64_t _table_size = 0;
    uint64_t _cache_size = 0;

    // Hash tables that hold data in memory, one per storer/consumer thread.
    sm_cache* _root_caches[MAX_STORERS];
    sm_stem_table* _stem_tables[MAX_STORERS];
    sm_root_table* _root_tables[MAX_STORERS];

    // Message queues between loader threads and storer threads. One SPSC
    // queue per loader/storer pair.
    sm_queue* _queues[MAX_STORERS][MAX_LOADERS];

    bool _enable_prune = false;
    const prune* _prune;

    input_queue* _input_queue;

    // Signal end of loader threads.
    std::atomic<bool> _done{false};

    // Track table conversion IDs.
    std::atomic<int> _convert{0};

    void load(int lid);
    void load_chunk(int lid, const sm_chunk &chunk);
    inline void load_sub(int lid, const char* sub, int len,
                         sm_read_kind kind, sm_bulk_msg* bulks);

    void incr(int sid);
    inline void incr_key(int sid, sm_key stem, sm_stem_offset off);

    void convert();
    void convert_table_mem(int sid);
    void convert_table_stream(int sid);

    void prefilter_table(int sid);

    void dump();
    void dump_table(int sid);
    void dump_table_stem(int sid);

    void restore();
    void restore_table(int sid);

    void export_csv();
    void export_csv_table(int sid);

    void annotate();
    void annotate_sub(const char* sub, int pos, int len, std::ofstream &ofs);

    void stats();
};

#endif
