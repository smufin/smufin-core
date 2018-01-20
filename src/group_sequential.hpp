/*
 * Copyright © 2015-2018 Barcelona Supercomputing Center (BSC)
 *
 * This file is part of SMUFIN Core. SMUFIN Core is released under the SMUFIN
 * Public License, and may not be used except in compliance with it. This file
 * is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; see the SMUFIN Public License for more details. You should have
 * received a copy of the SMUFIN Public License along with this file. If not,
 * see <https://github.com/smufin/smufin-core/blob/master/COPYING>.
 *
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2018
 */

#ifndef __SM_GROUP_SEQUENTIAL_H__
#define __SM_GROUP_SEQUENTIAL_H__

#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <google/sparse_hash_map>

#include "common.hpp"
#include "group.hpp"
#include "stage.hpp"

typedef struct sm_read_code {
    uint16_t len = 0;
    uint64_t seq[ENCODED_READ_LEN] = {0};
} sm_read_code;

typedef google::sparse_hash_map<std::string, sm_read_code> seq_table;
typedef google::sparse_hash_map<std::string, std::string> k2i_table;

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
// stored as sm_read_codes so as to minimize memory usage.
//
// After performing the initial sequential iterations, populate threads are
// spawned. At this point all required data is already indexed in memory, and
// each thread simply turns a subset of the data into groups and writes
// results to «group.PID-GID.json», where PID is the partition ID, and GID is
// the group thread ID.
class group_sequential : public stage
{
public:
    group_sequential(const sm_config &conf);
    void run();
    void stats();

private:
    uint64_t _leads_size = 0;

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

    void encode_read(std::string& str, sm_read_code& read);
    void decode_read(sm_read_code& read, std::string& str);

    void select_candidate(int gid, std::string& sid, std::string& seq,
                          std::string& dseq, std::vector<int>& pos, int dir);

    void populate(int gid);
    void populate_index(int gid, const std::string& lid,
                        const std::vector<std::string>& kmers, int kind,
                        kmer_count& keep, kmer_count& drop);
};

#endif
