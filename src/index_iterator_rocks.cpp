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
    path << this->_conf.output_path_filter << "/index-"
         << sm::types[this->_type] << "-" << sm::sets[this->_set] << "."
         << this->_pid << "-" << this->_iid << ".rdb";
    cout << "Prepare iterator: " << path.str() << endl;

    rdb_handle rdb;
    open_index_part_iter(this->_conf, this->_type, this->_set, this->_pid,
                         this->_iid, rdb);

    _it = rdb.db->NewIterator(rocksdb::ReadOptions());
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
