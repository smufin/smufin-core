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
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2018
 */

#ifndef __SM_ROCKSDB_H__
#define __SM_ROCKSDB_H__

#include <algorithm>
#include <sstream>
#include <string>

#include <rocksdb/db.h>
#include <rocksdb/env.h>
#include <rocksdb/merge_operator.h>

#include "common.hpp"

using namespace rocksdb;

typedef struct rdb_handle {
    rocksdb::DB* db;
    std::vector<rocksdb::ColumnFamilyHandle*> cfs;
} rdb_handle;

void set_options_type(const sm_config &conf,
                      rocksdb::ColumnFamilyOptions &options, sm_idx_type type);

void open_index(const sm_config &conf, sm_idx_type type,
                const std::string &path, const std::string &conf_file,
                rdb_handle &rdb);
void open_index_ro(const sm_config &conf, sm_idx_type type,
                   const std::string &path, rdb_handle &rdb);

void open_index_part_load(const sm_config &conf, sm_idx_type type,
                          sm_idx_set set, int pid, int iid,
                          rdb_handle &rdb);
void open_index_part_iter(const sm_config &conf, sm_idx_type type,
                          sm_idx_set set, int pid, int iid,
                          rdb_handle &rdb);

void open_index_full_load(const sm_config &conf, sm_idx_type type,
                          sm_idx_set set, rdb_handle &rdb);
void open_index_full_iter(const sm_config &conf, sm_idx_type type,
                          sm_idx_set set, rdb_handle &rdb);

void open_index_full_read(const sm_config &conf, sm_idx_type type,
                          sm_idx_set set, rdb_handle &rdb);

void open_groups_part(const sm_config &conf, int gid, rdb_handle &rdb);
void open_groups(const sm_config &conf, const std::string &path,
                 const std::string &conf_file, rdb_handle &rdb);

void encode_pos(const sm_pos_bitmap &p, std::string &s);
sm_pos_bitmap decode_pos(const std::string &s);

class PositionsMapOperator : public AssociativeMergeOperator
{
public:
    virtual bool Merge(const Slice& key, const Slice* existing_value,
                       const Slice& value, std::string* new_value,
                       Logger* logger) const override
    {
        sm_pos_bitmap existing;
        sm_pos_bitmap oper;
        sm_pos_bitmap result;
        oper = decode_pos(value.ToString());
        if (existing_value) {
            existing = decode_pos(existing_value->ToString());
        }
        for (int i = 0; i < POS_LEN; i++) {
            result.a[i] = existing.a[i] | oper.a[i];
            result.b[i] = existing.b[i] | oper.b[i];
        }
        encode_pos(result, *new_value);
        return true;
    }

    virtual const char* Name() const override
    {
        return "PositionsMapOperator";
    }
};

class IDListOperator : public AssociativeMergeOperator
{
public:
    IDListOperator(const sm_config &conf) : _conf(conf) {};

    virtual bool Merge(const Slice& key, const Slice* existing_value,
                       const Slice& value, std::string* new_value,
                       Logger* logger) const override
    {
        std::string existing;
        std::string oper;
        oper = std::string(value.ToString().c_str());
        if (existing_value) {
            existing = std::string(existing_value->ToString().c_str());
        }

        std::stringstream s;

        // Check if this merge operation exceeds the maximum number of reads
        // per kmer. In order to avoid counting the reads for every single
        // operation, the length of the value string is measured first as a
        // fast heuristic. The exact number is only counted when the length of
        // the value is 10 times bigger than the maximum number of reads
        // (assuming read IDs are usually at least 10 characters long).
        if (existing.size() > _conf.max_filter_reads * 10) {
            int count = std::count(existing.begin(), existing.end(), ' ');
            if (count > _conf.max_filter_reads) {
                s << existing;
                *new_value = s.str();
                return true;
            }
        }

        s << existing << " " << oper;
        *new_value = s.str();
        return true;
    }

    virtual const char* Name() const override
    {
        return "IDListOperator";
    }

private:
    const sm_config &_conf;
};

#endif
