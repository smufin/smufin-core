#ifndef __SM_ROCKSDB_H__
#define __SM_ROCKSDB_H__

#include <sstream>
#include <string>

#include <rocksdb/db.h>
#include <rocksdb/env.h>
#include <rocksdb/merge_operator.h>

#include "common.hpp"

using namespace rocksdb;

void set_options_type(rocksdb::Options &options, std::string type);
void set_options_filter(rocksdb::Options &options);
void set_options_merge(rocksdb::Options &options);

void open_filter(rocksdb::DB** db, const sm_config &conf, std::string type,
                 sm_set set, int pid, bool ro = false);
void open_merge(rocksdb::DB** db, const sm_config &conf, std::string type,
                sm_set set, bool ro = false);

void encode_pos(const sm_pos_bitmap &p, std::string &s);
sm_pos_bitmap decode_pos(char const *s);

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
        oper = decode_pos(value.data());
        if (existing_value) {
            existing = decode_pos(existing_value->data());
        }
        result.a[0] = existing.a[0] | oper.a[0];
        result.a[1] = existing.a[1] | oper.a[1];
        result.b[0] = existing.b[0] | oper.b[0];
        result.b[1] = existing.b[1] | oper.b[1];
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
    virtual bool Merge(const Slice& key, const Slice* existing_value,
                       const Slice& value, std::string* new_value,
                       Logger* logger) const override
    {
        std::string existing;
        std::string oper;
        oper = std::string(value.data());
        if (existing_value) {
            existing = std::string(existing_value->data());
        }
        std::stringstream s;
        s << existing << oper;
        *new_value = s.str();
        return true;
    }

    virtual const char* Name() const override
    {
        return "IDListOperator";
    }
};

#endif
