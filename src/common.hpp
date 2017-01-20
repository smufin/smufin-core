#ifndef __SM_COMMON_H__
#define __SM_COMMON_H__

#include <map>
#include <string>
#include <set>

#include <zlib.h>

#include "config.hpp"
#include "hash.hpp"
#include "kseq.h"

#define BASE_LEN 4
#define MAP_LEN 5
#define MAP_FILE_LEN 1024 // (BASE_LEN ^ MAP_LEN)

// Convert a string of 5 chars/bytes of the 4-base ACGT alphabet (40 bits)
// into a unique unsigned int identifier in the [0-1024) range (10 bits). The
// idea is to take the 2nd and 3rd least significant bits of each byte as
// follows:
//   A -> 65 -> 01000001 -> 00
//   C -> 67 -> 01000011 -> 01
//   G -> 71 -> 01000111 -> 11
//   T -> 84 -> 01010100 -> 10
#define hash_5mer(h) ({                   \
        (h) = ((h) & 25870861830) >> 1;     \
        (h) = ((h) & 0xFFFF) + ((h) >> 14); \
        (h) = ((h) & 0xFF) + ((h) >> 4);    \
        (h) = ((h) & 0xFF) + (((h) & 0xFF00) >> 6); })

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

#define NUM_TYPES 3
enum sm_idx_type {
    SEQ, // Sequence ID to sequence
    K2I, // Kmer to sequence IDs
    I2P, // ID to positions
};

#define NUM_SETS 3
enum sm_idx_set {
    NN, // Normal Non-mutated reads.
    TN, // Tumor Non-mutated reads.
    TM, // Tumor Mutated reads.
};

namespace sm {
    // Nucleotide alphabet, sorted and indexed by code.
    const char alpha[] = {'A', 'C', 'G', 'T'};

    // Map ASCII position to internal code. Note that partition maps use a
    // different encoding based on binary representation instead, see
    // hash_5mer.
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

    // Index and set names, sorted as the sm_index and sm_set enums,
    // respectively.
    const std::array<std::string, NUM_TYPES> types = {"seq", "k2i", "i2p"};
    const std::array<std::string, NUM_SETS> sets = {"nn", "tn", "tm"};

    // Map of filter indexes to valid set names.
    const std::map<sm_idx_type, std::set<sm_idx_set>> indexes = {
        {SEQ, {NN, TN, TM}},
        {K2I, {NN, TN, TM}},
        {I2P, {TM}}
    };
}

// Arrays that map which prefixes are to be processed on the current
// process (l1), and storer/consumer threads (l2), distributing them as
// evenly as possible according to MAP_FILE.
extern int map_l1[MAP_FILE_LEN];
extern int map_l2[MAP_FILE_LEN];

KSEQ_INIT(gzFile, gzread);

uint64_t strtob4(const char *str);
int lq_count(const char *str, int len);
void revcomp(char seq[], int len);

#endif
