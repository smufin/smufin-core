#include "filter_iterator_rocks.hpp"

#include <iostream>
#include <string>
#include <sstream>

#include "db.hpp"

using std::cout;
using std::endl;
using std::string;

template <typename T>
bool rocks_iterator<T>::init()
{
    std::ostringstream path;
    path << this->_conf.output_path << "/filter-" << sm::types[this->_type]
         << "-" << sm::sets[this->_set] << "." << this->_pid << ".rdb";
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

bool k2i_rocks_iterator::next()
{
    if (_it->Valid()) {
        _elem = new k2i_t(_it->key().ToString(), _it->value().ToString());
        _it->Next();
        return true;
    }
    return false;
}

bool i2p_rocks_iterator::next()
{
    if (_it->Valid()) {
        sm_pos_bitmap p = decode_pos(_it->value().data());
        _elem = new i2p_t(_it->key().ToString(), p);
        _it->Next();
        return true;
    }
    return false;
}
