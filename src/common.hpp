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

#ifndef __SM_COMMON_H__
#define __SM_COMMON_H__

#include <map>
#include <string>
#include <set>

#include <htslib/kseq.h>
#include <zlib.h>

#include "config.hpp"
#include "hash.hpp"

#define VERSION "2.0.0-b2"

#define CEIL(x,y) (((x) + (y) - 1) / (y))

#define BASE_LEN 4
#define MAP_LEN 6
#define MAP_FILE_LEN 4096 // (BASE_LEN ^ MAP_LEN)

#ifndef MAX_READ_LEN
#define MAX_READ_LEN 120
#endif

#define POS_LEN CEIL(MAX_READ_LEN, 64)

#define MAX_STORERS 128
#define MAX_LOADERS 128

// Convert a string of 6 chars/bytes of the 4-base ACGT alphabet (48 bits)
// into a unique unsigned int identifier in the [0-4096) range (12 bits). The
// idea is to take the 2nd and 3rd least significant bits of each byte as
// follows:
//   A -> 65 -> 01000001 -> 00
//   C -> 67 -> 01000011 -> 01
//   G -> 71 -> 01000111 -> 11
//   T -> 84 -> 01010100 -> 10
#define map_mer(h) ({                   \
    (h) = ((h) & 6622940628486) >> 1;     \
    (h) = ((h) & 0xFFFFFF) + ((h) >> 22); \
    (h) = ((h) & 0xF) + (((h) & 0xF00) >> 4) + (((h) & 0xF0000) >> 8); });

typedef uint64_t sm_key;

enum sm_read_kind : uint8_t {
    NORMAL_READ, CANCER_READ
};

// Map positions of candidate kmers in a sequence, in both directions: A and
// B. There are POS_LEN 64-bit elements in each array/direction, allowing
// sequences of up to MAX_READ_LEN+k-1 bases. The first element in the array
// maps positions 0..63, while the second element maps 64..127, etc. E.g. The
// following sequence of 30 bases has two candidate 12-mers starting at
// positions 0 and 2 in direction A, and one candidate starting at position 3
// in direction B:
//   A: GGGGTGCAGGTCCAAGGAAAGTCTTAGTGT
//      GGGGTGCAGGTC (0)
//        GGTGCAGGTCCA (2)
//   B: ACACTAAGACTTTCCTTGGACCTGCACCCC
//         CTAAGACTTTCC (3)
// Which results in the following bitmap:
//   a[0]: 0000000000000000000000000000000000000000000000000000000000000101 = 5
//   b[0]: 0000000000000000000000000000000000000000000000000000000000001000 = 8
typedef struct sm_pos_bitmap {
    uint64_t a[POS_LEN] = {0};
    uint64_t b[POS_LEN] = {0};
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

enum sm_dir {
    DIR_A,
    DIR_B
};

namespace sm {
    // Nucleotide alphabet, sorted and indexed by code.
    const char alpha[] = {'A', 'C', 'G', 'T'};

    // Map ASCII position to internal code. Note that partition maps use a
    // different encoding based on binary representation instead, see
    // map_mer.
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

    // Encoded complementary mapping.
    const int comp_code[] = {3, 2, 1, 0};

    // Index and set names, sorted as the sm_index and sm_set enums,
    // respectively.
    const std::array<std::string, NUM_TYPES> types = {"seq", "k2i", "i2p"};
    const std::array<std::string, NUM_SETS> sets = {"nn", "tn", "tm"};

    // Map of filter indexes to valid set names.
    const std::map<sm_idx_type, std::set<sm_idx_set>> indexes = {
        {SEQ, {NN, TN, TM}},
        {K2I, {NN, TN}},
        {I2P, {TM}}
    };

    // Supported filter format names.
    const std::set<std::string> formats = {"plain", "rocks"};
}

// Arrays that map which prefixes are to be processed on the current
// process (l1), and storer/consumer threads (l2), distributing them as
// evenly as possible according to MAP_FILE.
extern int map_l1[MAP_FILE_LEN];
extern int map_l2[MAP_FILE_LEN];

KSEQ_INIT(gzFile, gzread);

int lq_count(const char *str, int len);

uint64_t strtob4(const char *str);
void b4tostr(uint64_t code, int len, char *str);

void rev(char seq[], int len);
void revcomp(char seq[], int len);
sm_key revcomp_code(sm_key key, int len);

sm_key to_root(sm_key stem, int len);
int min_order(char seq[], int len);

#endif
