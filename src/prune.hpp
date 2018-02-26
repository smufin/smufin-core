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

#ifndef __SM_PRUNE_H__
#define __SM_PRUNE_H__

#include <string>

#include <bf/all.hpp>
#include <concurrentqueue.h>
#include <readerwriterqueue.h>

#include "common.hpp"
#include "input.hpp"
#include "stage.hpp"

#define BULK_KEY_LEN 512
#define PRUNE_QUEUE_LEN 128

typedef struct {
    uint16_t num = 0;
    sm_key array[BULK_KEY_LEN];
} sm_bulk_key;

typedef moodycamel::ReaderWriterQueue<sm_bulk_key> sm_prune_queue;

class prune : public stage
{
public:
    prune(const sm_config &conf);
    void run();
    void stats();

    inline const bf::basic_bloom_filter* operator[](int sid) const {
        return _allowed[sid];
    };

private:
    uint64_t _all_size = 0;
    uint64_t _allowed_size = 0;

    bf::basic_bloom_filter* _all[MAX_STORERS];
    bf::basic_bloom_filter* _allowed[MAX_STORERS];

    sm_prune_queue* _queues[MAX_STORERS][MAX_LOADERS];

    input_queue* _input_queue;

    // Signal end of loader threads.
    std::atomic<bool> _done{false};

    void load(int lid);
    void load_chunk(int lid, const sm_chunk &chunk);
    inline void load_sub(int lid, const char* sub, int len,
                         sm_bulk_key* bulk);

    void add(int sid);
    void add_key(int sid, sm_key key);
};

#endif
