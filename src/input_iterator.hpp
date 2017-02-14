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
