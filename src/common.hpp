#ifndef __SM_COMMON_H__
#define __SM_COMMON_H__

#include <string>
#include <set>
#include <map>

#include <zlib.h>

#include <concurrentqueue.h>
#include <google/sparse_hash_map>

#include "config.hpp"
#include "hash.hpp"
#include "kseq.h"

#define MAX_STORERS 128
#define MAX_LOADERS 128

#define BASE_LEN 4
#define MAP_LEN 5
#define MAP_FILE_LEN 1024 // (BASE_LEN ^ MAP_LEN)

#define BULK_LEN 128
#define QMSG_LEN 512

// Convert a string of 5 chars/bytes of the 4-base ACGT alphabet (32 bits)
// into a unique unsigned int identifier in the [0-256) range (8 bits). The
// idea is to take the 2nd and 3rd least significant bits of each byte as
// follows:
//   A -> 65 -> 01000001 -> 00
//   C -> 67 -> 01000011 -> 01
//   G -> 71 -> 01000111 -> 11
//   T -> 84 -> 01010100 -> 10
#define hash_5c_map(h) ({                   \
        (h) = ((h) & 25870861830) >> 1;     \
        (h) = ((h) & 0xFFFF) + ((h) >> 14); \
        (h) = ((h) & 0xFF) + ((h) >> 4);    \
        (h) = ((h) & 0xFF) + (((h) & 0xFF00) >> 6); })

typedef uint64_t sm_key;

typedef struct sm_value {
    uint16_t v[4][4][2] = {{{0}}};
} sm_value;

typedef google::sparse_hash_map<sm_key, sm_value, sm_hasher<sm_key>> sm_table;
typedef google::sparse_hash_map<sm_key, uint8_t, sm_hasher<sm_key>> sm_cache;

enum sm_read_kind : uint8_t {
    NORMAL_READ, CANCER_READ
};

// Map positions of candidate kmers in a sequence, in both directions: A and
// B. There are two 64-bit elements in each array/direction, allowing
// sequences of up to 128+k-1 bases. The first element in the array maps
// positions 0..63, while the second element maps 64..127.  Note that a[1] and
// b[1] are always zero for seequences that don't require indexing more than
// 64 bases. E.g. The following sequence of 30 bases has two candidate
// 12-mers starting at positions 0 and 2 in direction A, and one candidate
// starting at position 3 in direction B:
//   A: GGGGTGCAGGTCCAAGGAAAGTCTTAGTGT
//      GGGGTGCAGGTC (0)
//        GGTGCAGGTCCA (2)
//   B: AGGGTGCAGGTCCAAGGAAAGTCTTAGTGT
//         GTGCAGGTCCAA (3)
// Which results in the following bitmap:
//   a[0]: 0000000000000000000000000000000000000000000000000000000000000101 = 5
//   a[1]: 0000000000000000000000000000000000000000000000000000000000000000 = 0
//   b[0]: 0000000000000000000000000000000000000000000000000000000000001000 = 8
//   b[1]: 0000000000000000000000000000000000000000000000000000000000000000 = 0
typedef struct sm_pos_bitmap {
    uint64_t a[2] = {0};
    uint64_t b[2] = {0};
} sm_pos_bitmap;

#define NUM_SETS 3
enum sm_set {
    NN, // Normal Non-mutated reads.
    TN, // Tumor Non-mutated reads.
    TM, // Tumor Mutated reads.
};

// Arrays that map which prefixes are to be processed on the current
// process (l1), and storer/consumer threads (l2), distributing them as
// evenly as possible according to MAP_FILE.
extern int map_l1[MAP_FILE_LEN];
extern int map_l2[MAP_FILE_LEN];

namespace sm {
    // Nucleotide alphabet, sorted and indexed by code.
    const char alpha[] = {'A', 'C', 'G', 'T'};

    // Map ASCII position to internal code. Note that partition maps use a
    // different encoding based on binary representation instead, see
    // hash_5c_map.
    //   A -> 0 (pos. 65)
    //   C -> 1 (pos. 67)
    //   G -> 2 (pos. 71)
    //   T -> 3 (pos. 84)
    const char code[] = "-------------------------------------------"
                        "----------------------0-1---2------------3";

    // Map ASCII position to complementary nucleotide.
    //   A -> T (pos. 65)
    //   C -> G (pos. 67)
    //   G -> C (pos. 71)
    //   N -> N (pos. 78)
    //   T -> A (pos. 84)
    const char comp[] = "-------------------------------------------"
                        "----------------------T-G---C------N-----A";

    // Set names, sorted as the sm_set enum.
    const std::array<std::string, NUM_SETS> sets = {"nn", "tn", "tm"};

    // Map of filter indexes to valid set names.
    const std::map<std::string, std::set<std::string>> types = {
        {"seq", {"nn", "tn", "tm"}},
        {"k2i", {"nn", "tn", "tm"}},
        {"i2p", {"tm"}}
    };
}

KSEQ_INIT(gzFile, gzread);

uint64_t strtob4(const char *str);
int lq_count(const char *str, int len);
void revcomp(char seq[], int len);
void spawn(std::string name, std::function<void(int)> func, int n);
float estimate_sparse(uint64_t n, size_t k, size_t v);

#endif
