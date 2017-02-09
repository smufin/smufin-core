#ifndef __SM_INPUT_ITERATOR_FASTQ_H__
#define __SM_INPUT_ITERATOR_FASTQ_H__

#include "common.hpp"

#include <concurrentqueue.h>

#define MAX_SPLITS 10

typedef struct {
    kseq_t *seq;
    sm_read_kind kind;
    int num_splits;
    int splits[MAX_SPLITS][2] = {{0}};
} sm_split_read;

class input_queue
{
public:
    input_queue(const sm_config &conf);
    bool try_dequeue(std::string &file);

    std::atomic<int> len{0};

private:
    const sm_config &_conf;

    // SPMC queue to be initialized at startup time with the list of input
    // files to be processed. Idle producer threads will can try to read from
    // the queue until there are no input files left.
    moodycamel::ConcurrentQueue<std::string> _queue;
};

class input_iterator_fastq
{
public:
    input_iterator_fastq(const sm_config &conf) : _conf(conf) {};
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
