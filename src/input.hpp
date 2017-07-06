#ifndef __SM_INPUT__H__
#define __SM_INPUT__H__

#include <functional>

#include <concurrentqueue.h>

#include "common.hpp"

// Maximum number of allowed splits in a read. This works for the current
// range of read and kmer lengths, but needs to be increased for longer reads
// or shorter kmers.
#define MAX_SPLITS 10

// Internal read representation used for iteration over input reads; in
// addition to unique ID, sequence & qualities, and length of the sequence, it
// can also include splits to identify sub-sequences separated by undefined
// bases, which are handled as different in the smufin pipeline.
typedef struct {
    char *id;
    char *seq;
    char *qual;
    int len;
    int num_splits;
    int splits[MAX_SPLITS][2] = {{0}};
} sm_read;

// A chunk of the input to be processed at a time by a loader thread.
typedef struct {
    std::string file;
    uint64_t begin;
    uint64_t end;
    sm_read_kind kind;
} sm_chunk;

// Simple input queue that splits each file to be processed as a single chunk.
// Works for any kind of input format, but paralellization is constrained by
// number of input files.
class input_queue
{
public:
    input_queue(const sm_config &conf) : _conf(conf) {};
    void init();
    bool try_dequeue(sm_chunk &chunk);

    std::atomic<int> len{0};

    template<typename T> static input_queue* create(const sm_config &conf)
    {
        return new T(conf);
    };

protected:
    const sm_config &_conf;

    // SPMC queue to be initialized at startup time with the list of input
    // chunks to be processed. Idle producer threads will try to read from
    // the queue until there are no chunks left.
    moodycamel::ConcurrentQueue<sm_chunk> _queue;
};

// BAM-only input queue that splits every single file into «conf.num_loaders»
// chunks. The generation of chunks depends on BAI indexes to find valid BAM
// addresses, and may not be completely balanced for smaller BAM files.
class input_queue_bam_chunks : public input_queue
{
public:
    input_queue_bam_chunks(const sm_config &conf) : input_queue(conf) {};
    void init();

private:
    bool chunk_bam(const std::string bam_file, const int num_chunks,
                   std::vector<uint64_t> &offsets);
};

typedef std::function<input_queue*(const sm_config &conf)> input_queue_s;

#endif
