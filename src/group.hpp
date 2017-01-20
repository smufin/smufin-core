#ifndef __SM_GROUP_H__
#define __SM_GROUP_H__

#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <google/sparse_hash_map>

#include "common.hpp"
#include "stage.hpp"

#define RMAX 100
#define KMIN 0
#define KMAX 100

#define MAX_GROUPERS 128

#define ENCODED_READ_LEN 4 // ceil(max_read_len, 64/2)

typedef std::array<std::vector<int>, 2> p_value;
typedef std::array<std::vector<std::string>, 2> k_value;
typedef std::array<std::unordered_set<std::string>, 2> i_value;

// l2p: Lead ID to positions, direction A [0] and B [1]
// l2k: Lead ID to kmers, direction A [0] and B [1]
// l2i: Lead ID to sequence IDs, normal N [0] and tumoral T [1]
// l2r: Lead ID to lead sequence
typedef google::sparse_hash_map<std::string, p_value> l2p_table;
typedef google::sparse_hash_map<std::string, k_value> l2k_table;
typedef google::sparse_hash_map<std::string, i_value> l2i_table;
typedef google::sparse_hash_map<std::string, std::string> l2r_table;

typedef struct sm_read {
    uint16_t len = 0;
    uint64_t seq[ENCODED_READ_LEN] = {0};
} sm_read;

typedef google::sparse_hash_map<std::string, sm_read> seq_table;
typedef google::sparse_hash_map<std::string, std::string> k2i_table;

typedef std::array<std::unordered_map<std::string, int>, 2> kmer_count;

// Group stage initially designed to be able to run on MN3. There is a focus
// on not exceeded a certain amount of memory, and all reads from the merged
// RocksDB indexes are performed iterating sequentially so as to maximize
// performance retrieving data from a distributed GPFS filesystem.
//
// After loading shared data, work is split into different partitions/threads
// that work independently generating multiple files with the resulting
// groups. (Note that the number of partitions in the group stage can be
// different than in previous stages.)
//
// More specifically, there is an initial iteration over the i2p_tm index,
// which stores candidate leaders and positions, fetching also the leader's
// read from seq_tm, and ignoring those that aren't part of the current
// partition. Candidate leaders from i2p_tm that meet certain criteria are
// turned into effective group leaders and are initialized in l2r, l2p, and
// l2k; candidate kmers from the group's leader are also initialized as an
// empty string in _k2i. The second step involves iterating over the
// k2i_{nn,tn} indexes, populating the empty values in _k2i with real data
// from the RocksDB databases. And finally, there's another iteration over the
// seq_{nn,tn} indexes, loading sequences into _seq; sequences are encoded and
// stored as sm_reads so as to minimize memory usage.
//
// After performing the initial sequential iterations, populate threads are
// spawned. At this point all required data is already indexed in memory, and
// each thread simply turns a subset of the data into groups and writes
// results to «group.PID-GID.json», where PID is the partition ID, and GID is
// the group thread ID.
class group : public stage
{
public:
    group(const sm_config &conf);
    void run();
    void stats();

private:
    l2p_table* _l2p[MAX_GROUPERS];
    l2k_table* _l2k[MAX_GROUPERS];
    l2i_table* _l2i[MAX_GROUPERS];
    l2r_table* _l2r[MAX_GROUPERS];

    seq_table* _seq[2];
    k2i_table* _k2i[2];

    int _group_map_l1[MAP_FILE_LEN] = {0};
    int _group_map_l2[MAP_FILE_LEN] = {0};

    // Number of groups successfully generated by each grouper thread.
    uint64_t _num_groups[MAX_GROUPERS] = {0};

    void encode_read(std::string& str, sm_read& read);
    void decode_read(sm_read& read, std::string& str);
    void get_positions_a(uint64_t bitmap[2], std::vector<int> *pos);
    void get_positions_b(uint64_t bitmap[2], std::vector<int> *pos, int len);

    bool match_window(std::vector<int> pos);
    void select_candidate(int gid, std::string& sid, std::string& seq,
                          std::string& dseq, std::vector<int>& pos, int dir);

    void populate(int gid);
    void populate_index(int gid, std::string& lid,
                        const std::vector<std::string>& kmers, int kind,
                        kmer_count& keep, kmer_count& drop);
};

#endif
