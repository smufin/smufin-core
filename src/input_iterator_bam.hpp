#ifndef __SM_INPUT_ITERATOR_BAM_H__
#define __SM_INPUT_ITERATOR_BAM_H__

#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include "common.hpp"
#include "input.hpp"

class input_iterator_bam
{
public:
    input_iterator_bam(const sm_config &conf, const sm_chunk &chunk);
    bool next(sm_read *read);

private:
    const sm_config &_conf;
    const sm_chunk &_chunk;

    samFile *_in;
    bam1_t *_record;
    char _id[256];
    char _seq[256];
    char _qual[256];
};

#endif
