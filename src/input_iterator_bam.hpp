#ifndef __SM_INPUT_ITERATOR_BAM_H__
#define __SM_INPUT_ITERATOR_BAM_H__

#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include "common.hpp"
#include "input.hpp"
#include "input_iterator.hpp"

class input_iterator_bam : public input_iterator
{
public:
    input_iterator_bam(const sm_config &conf, const sm_chunk &chunk);
    bool next(sm_read *read);

private:
    samFile *_in;
    bam1_t *_record;
    char _id[256];
    char _seq[256];
    char _qual[256];
};

#endif
