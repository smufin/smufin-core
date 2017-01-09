#include "db.hpp"

#include <sstream>
#include <string>

#include "common.hpp"

rocksdb::Options get_rocks_options(std::string type)
{
    rocksdb::Options options;
    options.create_if_missing = true;

    options.disableDataSync = true;
    options.disable_auto_compactions = true;

    options.env->SetBackgroundThreads(4, Env::Priority::HIGH);
    options.env->SetBackgroundThreads(2, Env::Priority::LOW);

    options.num_levels = 2;
    options.level0_file_num_compaction_trigger = -1;
    options.level0_slowdown_writes_trigger = -1;
    options.level0_stop_writes_trigger = -1;
    options.soft_pending_compaction_bytes_limit = 0;
    options.hard_pending_compaction_bytes_limit = 0;
    options.write_buffer_size = 64UL * 1024 * 1024;
    options.max_write_buffer_number = 8;
    options.min_write_buffer_number_to_merge = 1;
    options.max_background_flushes = 4;
    options.max_background_compactions = 2;
    options.base_background_compactions = 2;
    options.target_file_size_base = 1024UL * 1024 * 1024;

    options.WAL_ttl_seconds = 0;
    options.WAL_size_limit_MB = 0;

    options.compression = rocksdb::kLZ4Compression;
    options.max_open_files = -1;

    if (type == "k2i") {
        options.merge_operator.reset(new IDListOperator());
    }

    if (type == "i2p") {
        options.merge_operator.reset(new PositionsMapOperator());
    }

    return options;
}

void encode_pos(const sm_pos_bitmap &p, std::string &s)
{
    std::stringstream e;
    e << std::hex << p.a[0] << " " << std::hex << p.a[1] << " "
      << std::hex << p.b[0] << " " << std::hex << p.b[1];
    s = e.str();
}

sm_pos_bitmap decode_pos(char const *s)
{
   sm_pos_bitmap p;
   std::istringstream(s)
       >> std::hex >> p.a[0] >> std::hex >> p.a[1]
       >> std::hex >> p.b[0] >> std::hex >> p.b[1];
   return p;
}
