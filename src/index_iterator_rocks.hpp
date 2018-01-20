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

#ifndef __SM_INDEX_ITERATOR_ROCKS_H__
#define __SM_INDEX_ITERATOR_ROCKS_H__

#include <rocksdb/db.h>

#include "common.hpp"
#include "index_iterator.hpp"

template <typename T>
class rocks_iterator : public index_iterator<T>
{
public:
    rocks_iterator(const sm_config &conf, sm_idx_set set, int pid, int iid,
                   sm_idx_type type)
        : index_iterator<T>(conf, set, pid, iid), _type(type) {};

    bool init();

protected:
    rocksdb::Iterator* _it;
    sm_idx_type _type;
};

class seq_rocks_iterator : public rocks_iterator<seq_t>
{
public:
    seq_rocks_iterator(const sm_config &conf, sm_idx_set set, int pid, int iid)
        : rocks_iterator<seq_t>(conf, set, pid, iid, SEQ) {};
    bool next();
};

class k2i_rocks_iterator : public rocks_iterator<k2i_t>
{
public:
    k2i_rocks_iterator(const sm_config &conf, sm_idx_set set, int pid, int iid)
        : rocks_iterator<k2i_t>(conf, set, pid, iid, K2I) {};
    bool next();
};

class i2p_rocks_iterator : public rocks_iterator<i2p_t>
{
public:
    i2p_rocks_iterator(const sm_config &conf, sm_idx_set set, int pid, int iid)
        : rocks_iterator<i2p_t>(conf, set, pid, iid, I2P) {};
    bool next();
};

#endif
