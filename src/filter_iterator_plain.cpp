#include "filter_iterator_plain.hpp"

#include <iostream>
#include <string>
#include <sstream>

using std::cout;
using std::endl;
using std::string;

template <typename T>
bool plain_iterator<T>::init()
{
    std::ostringstream file;
    file << this->_conf.output_path << "/filter-" << sm::types[this->_type]
         << "-" << sm::sets[this->_set] << "." << this->_pid << ".txt";
    cout << "Prepare iterator: " << file.str() << endl;
    this->_in.open(file.str());
    if (!this->_in.good()) {
        cout << "Failed to open: " << file.str() << endl;
        return false;
    }
    this->_elem = nullptr;
    return true;
}

bool seq_plain_iterator::next()
{
    string id;
    string seq;
    if (_in >> id >> seq) {
        delete _elem;
        _elem = new seq_t(id, seq);
        return true;
    }
    return false;
}

bool k2i_plain_iterator::next()
{
    string kmer;
    int len = 0;
    if (_in >> kmer >> len) {
        std::stringstream s;
        for (int i = 0; i < len; i++) {
            string sid;
            _in >> sid;
            s << sid << " ";
        }
        delete _elem;
        _elem = new k2i_t(kmer, s.str());
        return true;
    }
    return false;
}

bool i2p_plain_iterator::next()
{
    string id;
    sm_pos_bitmap p;
    if (_in >> id) {
        for (int i = 0; i < POS_LEN; i++)
            _in >> p.a[i];
        for (int i = 0; i < POS_LEN; i++)
            _in >> p.b[i];
        delete _elem;
        _elem = new i2p_t(id, p);
        return true;
    }
    return false;
}
