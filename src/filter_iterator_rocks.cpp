#include "filter_iterator_rocks.hpp"

#include <iostream>
#include <string>
#include <sstream>

using std::cout;
using std::endl;
using std::string;

bool seq_rocks_iterator::init()
{
    std::ostringstream path;
    path << _conf.output_path << "/filter-seq-" << _set << "." << _pid << ".rdb";
    cout << "Prepare iterator: " << path.str() << endl;

    rocksdb::DB* db;
    rocksdb::Options options;
    rocksdb::Status status;

    status = rocksdb::DB::OpenForReadOnly(options, path.str(), &db);
    if (!status.ok())
        return false;

    _it = db->NewIterator(rocksdb::ReadOptions());
    _it->SeekToFirst();
    return true;
}

bool seq_rocks_iterator::next()
{
    if (_it->Valid()) {
        _elem = new seq_t(_it->key().ToString(), _it->value().ToString());
        _it->Next();
        return true;
    }
    return false;
}
