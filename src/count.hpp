#ifndef __SM_COUNT_H__
#define __SM_COUNT_H__

#include <string>

#include <folly/ProducerConsumerQueue.h>
#include <google/sparse_hash_map>

#include "common.hpp"
#include "input.hpp"
#include "prune.hpp"
#include "stage.hpp"

#define BULK_MSG_LEN 128
#define COUNT_QUEUE_LEN 512

// A value of the hashtable that contains normal and tumoral counters for all
// inflections of a stem. The multidimensional array `v' is indexed as
// follows: v[A][B][C], where A and B are the numeric code of the first and
// last bases (as defined by sm::code), and C is the kind (as defined by
// sm_read_kind).
typedef struct sm_value {
    uint16_t v[4][4][2] = {{{0}}};
} sm_value;

typedef google::sparse_hash_map<sm_key, sm_value, sm_hasher<sm_key>> sm_table;
typedef google::sparse_hash_map<sm_key, uint8_t, sm_hasher<sm_key>> sm_cache;

// Contains data to calculate a position within a sm_value multidimensional
// array. `first' and `last' are integers in the range 0..3 and contain codes
// representing the first and last bases of a kmer using the same mapping as
// sm::code.
typedef struct {
    uint8_t first;
    uint8_t last;
    sm_read_kind kind;
} sm_value_offset;

typedef std::pair<sm_key, sm_value_offset> sm_msg;

typedef struct {
    uint16_t num = 0;
    sm_msg array[BULK_MSG_LEN];
} sm_bulk_msg;

typedef folly::ProducerConsumerQueue<sm_bulk_msg> sm_queue;

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

    inline const sm_table* operator[](int sid) const { return _tables[sid]; };

private:
    uint64_t _table_size = 0;
    uint64_t _cache_size = 0;

    // Hash tables that hold data in memory, one per storer/consumer thread.
    sm_table* _tables[MAX_STORERS];
    sm_cache* _caches[MAX_STORERS];

    // Message queues between loader threads and storer threads. One SPSC
    // queue per loader/storer pair.
    sm_queue* _queues[MAX_STORERS][MAX_LOADERS];

    bool _enable_prune = false;
    const prune* _prune;

    input_queue* _input_queue;

    // Signal end of loader threads.
    std::atomic<bool> _done{false};

    void load(int lid);
    void load_chunk(int lid, const sm_chunk &chunk);
    inline void load_sub(int lid, const char* sub, int len,
                         sm_read_kind kind, sm_bulk_msg* bulks);

    void incr(int sid);
    inline void incr_key(int sid, sm_key key, sm_value_offset off);

    void dump();
    void dump_table(int sid);

    void restore();
    void restore_table(int sid);

    void export_csv();
    void export_csv_table(int sid);

    void stats();
};

#endif
