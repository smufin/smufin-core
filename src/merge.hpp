#ifndef __SM_MERGE_H__
#define __SM_MERGE_H__

#include <rocksdb/db.h>

#include "common.hpp"
#include "stage.hpp"

class merge : public stage
{
public:
    merge(const sm_config &conf);
    void run();
    void stats();

private:
    void load(std::string type, std::string set);
    void load_seq(rocksdb::DB* db, std::string set, int i);
    void load_k2i(rocksdb::DB* db, std::string set, int i);
    void load_i2p(rocksdb::DB* db, std::string set, int i);
};

#endif
