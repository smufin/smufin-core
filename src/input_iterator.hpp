#ifndef __SM_INPUT_ITERATOR_H__
#define __SM_INPUT_ITERATOR_H__

#include "common.hpp"

#define MAX_SPLITS 10

typedef struct {
    kseq_t *seq;
    sm_read_kind kind;
    int num_splits;
    int splits[MAX_SPLITS][2] = {{0}};
} sm_split_read;

class input_iterator
{
public:
    input_iterator(const sm_config &conf) : _conf(conf) {};
    bool init(std::string file);
    bool next(sm_split_read *read);

private:
    const sm_config &_conf;

    std::string _file;
    sm_read_kind _kind;
    gzFile _in;
    kseq_t *_seq;
};

#endif
