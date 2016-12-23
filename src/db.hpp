#ifndef __SM_ROCKSDB_H__
#define __SM_ROCKSDB_H__

#include <sstream>
#include <string>

#include <rocksdb/db.h>
#include <rocksdb/merge_operator.h>

#include "common.hpp"

using namespace rocksdb;

rocksdb::Options get_rocks_options(std::string type);
void encode_pos(std::string &s, sm_pos_bitmap &p);
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
        encode_pos(*new_value, result);
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
