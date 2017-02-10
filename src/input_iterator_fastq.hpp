#ifndef __SM_INPUT_ITERATOR_FASTQ_H__
#define __SM_INPUT_ITERATOR_FASTQ_H__

#include "common.hpp"
#include "input.hpp"

class input_iterator_fastq
{
public:
    input_iterator_fastq(const sm_config &conf, const sm_chunk &chunk);
    bool next(sm_split_read *read);

private:
    const sm_config &_conf;
    const sm_chunk &_chunk;

    gzFile _in;
    kseq_t *_seq;
};

#endif
