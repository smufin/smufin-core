#ifndef __SM_COMMON_H__
#define __SM_COMMON_H__

#include <hash.hpp>

#include <mutex>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <concurrentqueue.h>
#include <boost/atomic.hpp>
#include <google/sparse_hash_map>
#include <folly/ProducerConsumerQueue.h>

// Expected number of keys per storer/thread.
// #define TABLE_LEN 250000000
#define TABLE_LEN 300000000
#define CACHE_LEN 3000000000

#define NUM_STORERS 16
#define MAX_LOADERS 16

#define BASE_LEN 4
#define MAP_LEN 4
#define MAP_FILE_LEN 256 // (BASE_LEN ^ MAP_LEN)

#ifndef READ_LEN
#define READ_LEN 100
#endif

#ifndef KMER_LEN
#define KMER_LEN 30
#endif

#define IMER_LEN (KMER_LEN - 2)

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

typedef struct sm_value {
    uint16_t v[4][4][2] = {{{0}}};
} sm_value;

typedef google::sparse_hash_map<sm_key, sm_value, sm_hasher<sm_key>> sm_table;
typedef google::sparse_hash_map<sm_key, uint8_t, sm_hasher<sm_key>> sm_cache;

enum sm_read_kind : uint8_t {
    NORMAL_READ, CANCER_READ
};

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

// Use ASCII codes to index base 4 values for ACGT. The array is
// equivalent to the following map, only slightly faster since it avoids
// hashing, etc.
// > std::unordered_map<char, char> code = {
// >     { 'A', '0' }, // pos. 65
// >     { 'C', '1' }, // pos. 67
// >     { 'G', '2' }, // pos. 71
// >     { 'T', '3' }, // pos. 84
// > };
const char code[] = "-------------------------------------------"
                    "----------------------0-1---2------------3";

// ASCII indexed bases to generate complementary strands.
// Equivalent to the following map:
// > std::unordered_map<char, char> comp = {
// >     { 'A', 'T' }, // pos. 65
// >     { 'C', 'G' }, // pos. 67
// >     { 'G', 'C' }, // pos. 71
// >     { 'N', 'N' }, // pos. 78
// >     { 'T', 'A' }, // pos. 84
// > };
const char comp[] = "-------------------------------------------"
                    "----------------------T-G---C------N-----A";

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
extern sm_cache caches[NUM_STORERS];

// Message queues between loader threads and storer threads. One SPSC queue
// per loader/storer pair.
extern folly::ProducerConsumerQueue<sm_bulk>* queues[NUM_STORERS][MAX_LOADERS];

#define NUM_SETS 3
enum sm_set {
    NN, // Normal Non-mutated reads.
    TM, // Tumor Mutated reads.
    TN, // Tumor Non-mutated reads.
};

extern std::vector<std::string> set_names;
extern std::mutex filter_mutex[NUM_SETS];
extern std::unordered_set<std::string> filter_reads[NUM_SETS];
extern std::unordered_map<std::string, std::pair<std::vector<uint8_t>, std::vector<uint8_t>>> filter_i2p[NUM_SETS];
extern std::unordered_map<std::string, std::unordered_set<std::string>> filter_k2i[NUM_SETS];

enum noshort_options {
    O_DISABLE_FILTER, O_DISABLE_STATS
};

uint64_t strtob4(const char *str);
int lq_count(const char *str);
void krevcomp(char s[]);

#endif
