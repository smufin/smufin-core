#include "filter_iterator_plain.hpp"

#include <iostream>
#include <string>
#include <sstream>

using std::cout;
using std::endl;
using std::string;

bool seq_plain_iterator::init()
{
    std::ostringstream file;
    file << _conf.output_path << "/filter-seq-" << _set << "." << _pid
         << ".txt";
    cout << "Prepare iterator: " << file.str() << endl;
    _in.open(file.str());
    if (!_in.good()) {
        cout << "Failed to open: " << file.str() << endl;
        return false;
    }
    return true;
}

bool seq_plain_iterator::next()
{
    string id;
    string seq;
    if (_in >> id >> seq) {
        _elem = new seq_t(id, seq);
        return true;
    }
    return false;
}

bool k2i_plain_iterator::init()
{
    std::ostringstream file;
    file << _conf.output_path << "/filter-k2i-" << _set << "." << _pid
         << ".txt";
    cout << "Prepare iterator: " << file.str() << endl;
    _in.open(file.str());
    if (!_in.good()) {
        cout << "Failed to open: " << file.str() << endl;
        return false;
    }
    return true;
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
        _elem = new k2i_t(kmer, s.str());
        return true;
    }
    return false;
}

bool i2p_plain_iterator::init()
{
    std::ostringstream file;
    file << _conf.output_path << "/filter-i2p-" << _set << "." << _pid
         << ".txt";
    cout << "Prepare iterator: " << file.str() << endl;
    _in.open(file.str());
    if (!_in.good()) {
        cout << "Failed to open: " << file.str() << endl;
        return false;
    }
    return true;
}

bool i2p_plain_iterator::next()
{
    string id;
    sm_pos_bitmap p;
    if (_in >> id >> p.a[0] >> p.a[1] >> p.b[0] >> p.b[1]) {
        _elem = new i2p_t(id, p);
        return true;
    }
    return false;
}
