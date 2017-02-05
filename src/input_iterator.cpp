#include "input_iterator.hpp"

#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

bool input_iterator::init(string file)
{
    gzFile _in = gzopen(file.c_str(), "rb");

    // Identify read kind from file name.
    _kind = NORMAL_READ;
    std::size_t found = file.find("_T_");
    if (found != std::string::npos) {
        _kind = CANCER_READ;
    }

    _seq = kseq_init(_in);
    return true;
}

bool input_iterator::next(sm_split_read *read)
{
    while (kseq_read(_seq) >= 0) {
        if (lq_count(_seq->qual.s, _seq->qual.l) > _seq->qual.l / 10)
            continue;

        read->seq = _seq;
        read->kind = _kind;
        read->num_splits = 0;

        int p = 0;
        int l = _seq->seq.l;
        int n = _seq->seq.l;
        char *ps;

        while ((ps = (char*) memchr(&_seq->seq.s[p], 'N', l - p)) != NULL) {
            n = ps - &_seq->seq.s[p];
            if (n > 0) {
                read->splits[read->num_splits][0] = p;
                read->splits[read->num_splits][1] = n;
                read->num_splits++;
                p += n;
            }
            p++;
        }

        n = l - p;
        if (n > 0) {
            read->splits[read->num_splits][0] = p;
            read->splits[read->num_splits][1] = n;
            read->num_splits++;
        }

        return true;
    }

    return false;
}
