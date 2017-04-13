#include "index_iterator_rocks.hpp"

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
    path << this->_conf.output_path << "/index-" << sm::types[this->_type]
         << "-" << sm::sets[this->_set] << "." << this->_pid << "-"
         << this->_iid << ".rdb";
    cout << "Prepare iterator: " << path.str() << endl;

    rocksdb::DB* db;
    rocksdb::Options options;
    rocksdb::Status status;

    set_options_type(options, this->_type);
    status = rocksdb::DB::OpenForReadOnly(options, path.str(), &db);
    if (!status.ok()) {
        cout << "Failed to open: " << path.str() << endl;
        return false;
    }

    _it = db->NewIterator(rocksdb::ReadOptions());
    _it->SeekToFirst();
    this->_elem = nullptr;
    return true;
}

bool seq_rocks_iterator::next()
{
    if (_it->Valid()) {
        delete _elem;
        _elem = new seq_t(_it->key().ToString(), _it->value().ToString());
        _it->Next();
        return true;
    }
    return false;
}

bool k2i_rocks_iterator::next()
{
    if (_it->Valid()) {
        delete _elem;
        _elem = new k2i_t(_it->key().ToString(), _it->value().ToString());
        _it->Next();
        return true;
    }
    return false;
}

bool i2p_rocks_iterator::next()
{
    if (_it->Valid()) {
        delete _elem;
        sm_pos_bitmap p = decode_pos(_it->value().ToString());
        _elem = new i2p_t(_it->key().ToString(), p);
        _it->Next();
        return true;
    }
    return false;
}
