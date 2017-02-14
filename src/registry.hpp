#ifndef __SM_REGISTRY_H__
#define __SM_REGISTRY_H__

#include "common.hpp"
#include "input.hpp"
#include "input_iterator.hpp"
#include "input_iterator_bam.hpp"
#include "input_iterator_fastq.hpp"

namespace sm
{
    const std::map<std::string, input_queue_s> input_queues = {
        {"fastq", &input_queue::create<input_queue>},
        {"bam", &input_queue::create<input_queue_bam_chunks>}
    };

    const std::map<std::string, input_iterator_s> input_iterators = {
        {"fastq", &input_iterator::create<input_iterator_fastq>},
        {"bam", &input_iterator::create<input_iterator_bam>}
    };
}

#endif
