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

#ifndef __SM_INPUT_ITERATOR_H__
#define __SM_INPUT_ITERATOR_H__

#include "common.hpp"
#include "input.hpp"

class input_iterator
{
public:
    input_iterator(const sm_config &conf, const sm_chunk &chunk)
        : _conf(conf), _chunk(chunk) {};
    virtual bool next(sm_read *read) = 0;

    template<typename T>
    static input_iterator* create(const sm_config &conf, const sm_chunk &chunk)
    {
        return new T(conf, chunk);
    };

protected:
    const sm_config &_conf;
    const sm_chunk &_chunk;
};

typedef std::function<input_iterator*(const sm_config &conf, const sm_chunk
        &chunk)> input_iterator_s;

#endif
