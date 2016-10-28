#ifndef __SM_GROUP_H__
#define __SM_GROUP_H__

#include <unordered_map>
#include <unordered_set>

#include <google/sparse_hash_map>
#include <rocksdb/db.h>

#include "common.hpp"

#define RMAX 100
#define KMIN 0
#define KMAX 100
#define DROP 500

#define MAX_GROUPERS 8

#define ENCODED_READ_LEN 4 // ceil(max_read_len, 64/2)

typedef std::array<std::vector<int>, 2> p_value;
typedef std::array<std::vector<std::string>, 2> k_value;
typedef std::array<std::unordered_set<std::string>, 2> i_value;

typedef google::sparse_hash_map<std::string, p_value> l2p_table;
typedef google::sparse_hash_map<std::string, k_value> l2k_table;
typedef google::sparse_hash_map<std::string, i_value> l2i_table;
typedef google::sparse_hash_map<std::string, std::string> l2r_table;

typedef struct sm_read {
    uint16_t len = 0;
    uint64_t seq[ENCODED_READ_LEN] = {0};
} sm_read;

typedef google::sparse_hash_map<std::string, sm_read> i2r_table;
typedef google::sparse_hash_map<std::string, std::string> kmer_table;

typedef std::array<std::unordered_map<std::string, int>, 2> index_count;

extern l2p_table* l2p[MAX_GROUPERS];
extern l2k_table* l2k[MAX_GROUPERS];
extern l2i_table* l2i[MAX_GROUPERS];
extern l2r_table* l2r[MAX_GROUPERS];

extern i2r_table* i2r[2];
extern kmer_table* k2i[2];

void encode_read(std::string& str, sm_read& read);
void decode_read(sm_read& read, std::string& str);
void rrevcomp(char read[], int len);
void get_positions_a(uint64_t bitmap[2], std::vector<int> *pos);
void get_positions_b(uint64_t bitmap[2], std::vector<int> *pos, int len);

bool match_window(std::vector<int> pos);
void select_candidate(int gid, std::string& sid, std::string& seq,
                      std::string& dseq, std::vector<int>& pos, int dir);

void populate(int pid, int gid);
void populate_index(int gid, std::string& lid,
                    const std::vector<std::string>& kmers, int kind,
                    index_count& keep, index_count& drop);

#endif
