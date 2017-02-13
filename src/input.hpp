#ifndef __SM_INPUT__H__
#define __SM_INPUT__H__

#include <concurrentqueue.h>

#include "common.hpp"

#define MAX_SPLITS 10

typedef struct {
    kseq_t *seq;
    int num_splits;
    int splits[MAX_SPLITS][2] = {{0}};
} sm_split_read;

typedef struct {
    std::string file;
    uint64_t begin;
    uint64_t end;
    sm_read_kind kind;
} sm_chunk;

class input_queue
{
public:
    input_queue(const sm_config &conf) : _conf(conf) {};
    void init();
    bool try_dequeue(sm_chunk &chunk);

    std::atomic<int> len{0};

protected:
    const sm_config &_conf;

    // SPMC queue to be initialized at startup time with the list of input
    // chunks to be processed. Idle producer threads will try to read from
    // the queue until there are no chunks left.
    moodycamel::ConcurrentQueue<sm_chunk> _queue;
};

class input_queue_bam_chunks : public input_queue
{
public:
    input_queue_bam_chunks(const sm_config &conf) : input_queue(conf) {};
    void init();

private:
    bool chunk_bam(const std::string bam_file, const int num_chunks,
                   std::vector<uint64_t> &offsets);
};

#endif
