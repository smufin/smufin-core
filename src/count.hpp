#ifndef __SM_COUNT_H__
#define __SM_COUNT_H__

#include <string>

#include "common.hpp"
#include "stage.hpp"

class count : public stage
{
public:
    count(const sm_config &conf);
    void run();

    inline const sm_table* operator[](int sid) const { return _tables[sid]; };

private:
    uint64_t _table_size = 0;
    uint64_t _cache_size = 0;

    // Hash tables that hold data in memory, one per storer/consumer thread.
    sm_table* _tables[MAX_STORERS];

    // Message queues between loader threads and storer threads. One SPSC
    // queue per loader/storer pair.
    sm_queue* _queues[MAX_STORERS][MAX_LOADERS];

    // Signal end of loader threads.
    std::atomic<bool> _done{false};

    // SPMC queue to be initialized at startup time with the list of input
    // files to be processed. Idle producer threads will read from the queue
    // until there are no input files left.
    moodycamel::ConcurrentQueue<std::string> _input_queue;
    std::atomic<int> _input_len{0};

    void load(int lid);
    void load_file(int lid, std::string file);
    inline void load_sub(int lid, const char* sub, int len,
                         sm_read_kind kind, sm_bulk* bulks);

    void incr(int sid);
    inline void incr_key(int sid, sm_cache* cache, sm_key key,
                         sm_value_offset off);

    void dump();
    void dump_table(int sid);

    void restore();
    void restore_table(int sid);

    void stats();
};

#endif
