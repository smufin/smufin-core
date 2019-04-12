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

#ifndef __SM_INPUT_ITERATOR_BAM_H__
#define __SM_INPUT_ITERATOR_BAM_H__

#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include "common.hpp"
#include "input.hpp"
#include "input_iterator.hpp"

// Reserve a big enough multiple of 64 to hold each one of the fields of the
// read. _seq and _qual need to be at least MAX_READ_LEN + 1; _id's length is
// always unknown, but it's usually significantly shorter.
#define BAM_BUF_LEN (((MAX_READ_LEN / 64) + 1) * 64)

class input_iterator_bam : public input_iterator
{
public:
    input_iterator_bam(const sm_config &conf, const sm_chunk &chunk);
    bool next(sm_read *read);

private:
    samFile *_in;
    bam_hdr_t *_header;
    bam1_t *_record;
    char _id[BAM_BUF_LEN];
    char _seq[BAM_BUF_LEN];
    char _qual[BAM_BUF_LEN];
};

#endif
