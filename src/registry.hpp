#ifndef __SM_REGISTRY_H__
#define __SM_REGISTRY_H__

#include "common.hpp"

#include "stage.hpp"
#include "prune.hpp"
#include "count.hpp"
#include "filter.hpp"
#include "merge.hpp"
#include "group_sequential.hpp"

#include "input.hpp"
#include "input_iterator.hpp"
#include "input_iterator_bam.hpp"
#include "input_iterator_fastq.hpp"

#include "index_format.hpp"
#include "index_format_plain.hpp"
#include "index_format_rocks.hpp"
#include "index_iterator.hpp"
#include "index_iterator_plain.hpp"
#include "index_iterator_rocks.hpp"

namespace sm
{
    const std::map<std::string, stage_s> stages = {
        {"prune", &stage::create<prune>},
        {"count", &stage::create<count>},
        {"filter", &stage::create<filter>},
        {"merge", &stage::create<merge>},
        {"group", &stage::create<group_sequential>}
    };

    const std::map<std::string, input_queue_s> input_queues = {
        {"fastq", &input_queue::create<input_queue>},
        {"bam", &input_queue::create<input_queue_bam_chunks>}
    };

    const std::map<std::string, input_iterator_s> input_iterators = {
        {"fastq", &input_iterator::create<input_iterator_fastq>},
        {"bam", &input_iterator::create<input_iterator_bam>}
    };

    const std::map<std::string, index_format_s> index_formats = {
        {"plain", &index_format::create<index_format_plain>},
        {"rocks", &index_format::create<index_format_rocks>}
    };

    const std::map<std::string, seq_iterator_s> seq_iterators = {
        {"plain", &create_index_iterator<seq_plain_iterator>},
        {"rocks", &create_index_iterator<seq_rocks_iterator>}
    };

    const std::map<std::string, k2i_iterator_s> k2i_iterators = {
        {"plain", &create_index_iterator<k2i_plain_iterator>},
        {"rocks", &create_index_iterator<k2i_rocks_iterator>}
    };

    const std::map<std::string, i2p_iterator_s> i2p_iterators = {
        {"plain", &create_index_iterator<i2p_plain_iterator>},
        {"rocks", &create_index_iterator<i2p_rocks_iterator>}
    };
}

#endif
