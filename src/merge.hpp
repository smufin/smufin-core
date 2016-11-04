#ifndef __SM_MERGE_H__
#define __SM_MERGE_H__

#include <map>
#include <set>

#include <rocksdb/db.h>

#include "common.hpp"
#include "stage.hpp"

class merge : public stage
{
public:
    merge(const sm_config &conf);
    void run();

private:
    const std::map<std::string, std::set<std::string>> _types;

    void load(std::string type, std::string set);
    void load_seq(rocksdb::DB* db, std::string set, int i);
    void load_k2i(rocksdb::DB* db, std::string set, int i);
    void load_i2p(rocksdb::DB* db, std::string set, int i);
};

#endif
