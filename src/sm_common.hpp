#ifndef __SM_COMMON_H__
#define __SM_COMMON_H__

#include <sm_hash.hpp>

#include <mutex>
#include <string>
#include <unordered_set>
#include <concurrentqueue.h>
#include <boost/atomic.hpp>
#include <sparsehash/sparse_hash_map>
#include <folly/ProducerConsumerQueue.h>

// Expected number of keys per storer/thread.
#define TABLE_LEN 250000000

#define NUM_STORERS 16
#define MAX_NUM_LOADERS 16

#define BASE_LEN 4
#define MAP_LEN 4
#define MAP_FILE_LEN 256 // (BASE_LEN ^ MAP_LEN)

#define KMER_LEN 30
#define ZBUF_LEN 2592
#define BULK_LEN 128
#define QMSG_LEN 512

// Convert a string of 4 chars/bytes of the 4-base ACGT alphabet (32 bits)
// into a unique unsigned int identifier in the [0-256) range (8 bits). The
// idea is to take the 2nd and 3rd least significant bits of each byte as
// follows:
//   A -> 65 -> 01000001 -> 00
//   C -> 67 -> 01000011 -> 01
//   G -> 71 -> 01000111 -> 11
//   T -> 84 -> 01010100 -> 10
#define hash_4c_map(h) ({                   \
        (h) = ((h) & 101058054) >> 1;       \
        (h) = ((h) & 0xFFFF) + ((h) >> 14); \
        (h) = ((h) & 0xFF) + ((h) >> 4); })

typedef uint64_t sm_key;
typedef std::pair<uint32_t, uint32_t> sm_value;
typedef google::sparse_hash_map<sm_key, sm_value, sm_hasher<sm_key>> sm_table;

enum sm_read_kind {
    NORMAL_READ, CANCER_READ
};

typedef std::pair<sm_key, sm_read_kind> sm_msg;
typedef struct {
    uint16_t num = 0;
    sm_msg array[BULK_LEN];
} sm_bulk;

// Use ASCII codes to index base 4 values for ACGT. The array is
// equivalent to the following map, only slightly faster since it avoids
// hashing, etc.
// > std::unordered_map<char, char> codes = {
// >     { 'A', '0' }, // pos. 65
// >     { 'C', '1' }, // pos. 67
// >     { 'G', '2' }, // pos. 71
// >     { 'T', '3' }, // pos. 84
// > };
const char codes[] = "-------------------------------------------"
                     "----------------------0-1---2------------3";

// SPMC queue to be initialized at startup time with the list of input
// files to be processed. Idle producer threads will read from the queue
// until there are no input files left.
extern moodycamel::ConcurrentQueue<std::string> input_queue;
extern boost::atomic_int input_count;

// Signal end of processing and filtering threads.
extern boost::atomic<bool> process_done;

// Arrays that map which prefixes are to be processed on the current
// process (l1), and storer/consumer threads (l2), distributing them as
// evenly as possible according to MAP_FILE.
extern int map_l1[MAP_FILE_LEN];
extern int map_l2[MAP_FILE_LEN];

// Hash tables that hold data in memory, one per storer/consumer thread.
extern sm_table tables[NUM_STORERS];

// Message queues between loader threads and storer threads. One SPSC queue
// per loader/storer pair.
extern folly::ProducerConsumerQueue<sm_bulk>* queues[NUM_STORERS][MAX_NUM_LOADERS];

extern std::unordered_set<std::string> filter_reads;
extern std::mutex filter_mutex;

enum noshort_options {
    O_DISABLE_FILTER, O_DISABLE_STATS
};

unsigned long int strtob4(const char *str);

inline char encode(char c)
{
    return codes[c];
}

inline int lq_count(char *str)
{
    int lq = 0;
    for (int i = 0; i < 80; i++) {
        int phred = str[i] - 33;
        if (phred < 20)
            lq++;
    }
    return lq;
}

#endif
