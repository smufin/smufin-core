#ifndef __SM_COUNT_H__
#define __SM_COUNT_H__

#include <string>

#include <folly/ProducerConsumerQueue.h>

#include "common.hpp"
#include "stage.hpp"

typedef uint64_t sm_key;

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
    sm_msg array[BULK_LEN];
} sm_bulk;

typedef folly::ProducerConsumerQueue<sm_bulk> sm_queue;

// Stage that reads input files, splits sequences into kmers, and builds a
// table of normal and tumoral kmer frequencies. `count' provides an in-memory
// implementation, and uses a cache that holds kmers that are seen only once.
// Frequency tables are indexed by stem instead of kmers, and each entry
// contains normal and tumoral counters for all inflections.
class count : public stage
{
public:
    count(const sm_config &conf);
    void run();

    inline const sm_table* operator[](int sid) const { return _tables[sid]; };

private:
    uint64_t _table_size = 0;
    uint64_t _cache_size = 0;

    // Hash tables that hold data in memory, one per storer/consumer thread.
    sm_table* _tables[MAX_STORERS];

    // Message queues between loader threads and storer threads. One SPSC
    // queue per loader/storer pair.
    sm_queue* _queues[MAX_STORERS][MAX_LOADERS];

    // Signal end of loader threads.
    std::atomic<bool> _done{false};

    // SPMC queue to be initialized at startup time with the list of input
    // files to be processed. Idle producer threads will read from the queue
    // until there are no input files left.
    moodycamel::ConcurrentQueue<std::string> _input_queue;
    std::atomic<int> _input_len{0};

    void load(int lid);
    void load_file(int lid, std::string file);
    inline void load_sub(int lid, const char* sub, int len,
                         sm_read_kind kind, sm_bulk* bulks);

    void incr(int sid);
    inline void incr_key(int sid, sm_cache* cache, sm_key key,
                         sm_value_offset off);

    void dump();
    void dump_table(int sid);

    void restore();
    void restore_table(int sid);

    void stats();
};

#endif
