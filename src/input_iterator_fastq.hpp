#ifndef __SM_INPUT_ITERATOR_FASTQ_H__
#define __SM_INPUT_ITERATOR_FASTQ_H__

#include "common.hpp"
#include "input.hpp"
#include "input_iterator.hpp"

class input_iterator_fastq : public input_iterator
{
public:
    input_iterator_fastq(const sm_config &conf, const sm_chunk &chunk);
    bool next(sm_read *read);
    bool _check = true;

private:
    gzFile _in;
    kseq_t *_seq;
};

#endif
