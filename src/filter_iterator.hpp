#ifndef __SM_FILTER_ITERATOR_H__
#define __SM_FILTER_ITERATOR_H__

#include <string>

#include "filter_format.hpp"

typedef std::pair<std::string, std::string> seq_t;
typedef std::pair<std::string, std::string> k2i_t;
typedef std::pair<std::string, sm_pos_bitmap> i2p_t;

template<typename T>
class filter_iterator
{
public:
    filter_iterator(const sm_config &conf, sm_idx_set set, int pid, int fid)
        : _conf(conf), _set(set), _pid(pid), _fid(fid) {};

    virtual bool init() = 0;
    virtual bool next() = 0;
    const T* get() { return _elem; };

protected:
    const sm_config &_conf;
    const sm_idx_set _set;
    const int _pid;
    const int _fid;
    const T* _elem;
};

template<typename T>
T* create_filter_iterator(const sm_config &conf, sm_idx_set set, int pid,
                          int fid)
{
    return new T(conf, set, pid, fid);
};

typedef std::function<filter_iterator<seq_t>*(const sm_config &conf,
        sm_idx_set set, int pid, int fid)> seq_iterator_s;
typedef std::function<filter_iterator<k2i_t>*(const sm_config &conf,
        sm_idx_set set, int pid, int fid)> k2i_iterator_s;
typedef std::function<filter_iterator<i2p_t>*(const sm_config &conf,
        sm_idx_set set, int pid, int fid)> i2p_iterator_s;

#endif
