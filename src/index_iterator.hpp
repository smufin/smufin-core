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

#ifndef __SM_INDEX_ITERATOR_H__
#define __SM_INDEX_ITERATOR_H__

#include <string>

#include "index_format.hpp"

typedef std::pair<std::string, std::string> seq_t;
typedef std::pair<std::string, std::string> k2i_t;
typedef std::pair<std::string, sm_pos_bitmap> i2p_t;

template<typename T>
class index_iterator
{
public:
    index_iterator(const sm_config &conf, sm_idx_set set, int pid, int iid)
        : _conf(conf), _set(set), _pid(pid), _iid(iid) {};

    virtual bool init() = 0;
    virtual bool next() = 0;
    const T* get() { return _elem; };

protected:
    const sm_config &_conf;
    const sm_idx_set _set;
    const int _pid;
    const int _iid;
    const T* _elem;
};

template<typename T>
T* create_index_iterator(const sm_config &conf, sm_idx_set set, int pid,
                          int iid)
{
    return new T(conf, set, pid, iid);
};

typedef std::function<index_iterator<seq_t>*(const sm_config &conf,
        sm_idx_set set, int pid, int iid)> seq_iterator_s;
typedef std::function<index_iterator<k2i_t>*(const sm_config &conf,
        sm_idx_set set, int pid, int iid)> k2i_iterator_s;
typedef std::function<index_iterator<i2p_t>*(const sm_config &conf,
        sm_idx_set set, int pid, int iid)> i2p_iterator_s;

#endif
