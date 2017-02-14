#include "input_iterator_fastq.hpp"

#include <fstream>
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

input_iterator_fastq::input_iterator_fastq(const sm_config &conf,
                                           const sm_chunk &chunk)
    : input_iterator(conf, chunk)
{
    gzFile _in = gzopen(_chunk.file.c_str(), "rb");
    _seq = kseq_init(_in);
}

bool input_iterator_fastq::next(sm_read *read)
{
    while (kseq_read(_seq) >= 0) {
        if (lq_count(_seq->qual.s, _seq->qual.l) > _seq->qual.l / 10)
            continue;

        read->id = _seq->name.s;
        read->seq = _seq->seq.s;
        read->qual = _seq->qual.s;
        read->len = _seq->seq.l;
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
