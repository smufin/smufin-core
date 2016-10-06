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

typedef std::array<std::vector<int>, 2> pos_value;
typedef std::array<std::vector<std::string>, 2> kmer_value;
typedef std::array<std::unordered_set<std::string>, 2> id_value;

typedef google::sparse_hash_map<std::string, pos_value> l2p_table;
typedef google::sparse_hash_map<std::string, kmer_value> l2k_table;
typedef google::sparse_hash_map<std::string, id_value> l2i_table;

typedef std::array<std::unordered_map<std::string, int>, 2> index_count;

extern l2p_table* l2p[MAX_GROUPERS];
extern l2k_table* l2k[MAX_GROUPERS];
extern l2i_table* l2i[MAX_GROUPERS];

extern rocksdb::DB* i2r[NUM_SETS];
extern rocksdb::DB* k2i[NUM_SETS];

void rrevcomp(char read[], int len);
void get_positions_a(uint64_t bitmap[2], std::vector<int> *pos);
void get_positions_b(uint64_t bitmap[2], std::vector<int> *pos, int len);

bool match_window(std::vector<int> pos);
void select_candidate(int gid, std::string sid, std::string seq,
                      std::vector<int>& pos, int dir);

void populate(int pid, int gid);
void populate_index(int gid, std::string& lid,
                    const std::vector<std::string>& kmers, int kind,
                    index_count& keep, index_count& drop);

#endif
