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
    file << _conf.output_path << "/filter-seq-" << _set << "." << _pid << ".txt";
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
