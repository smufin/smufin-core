#ifndef __SM_MERGE_H__
#define __SM_MERGE_H__

#include <rocksdb/db.h>

void load_i2p(rocksdb::DB* db, std::string path, std::string set, int i);
void load_k2i(rocksdb::DB* db, std::string path, std::string set, int i);
void load_seq(rocksdb::DB* db, std::string path, std::string set, int i);

#endif
