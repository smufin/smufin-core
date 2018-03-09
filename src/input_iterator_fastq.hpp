/*
 * Copyright © 2015-2018 Barcelona Supercomputing Center (BSC)
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

private:
    gzFile _in;
    kseq_t *_seq;
};

#endif
