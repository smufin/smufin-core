#include "db.hpp"

#include <sstream>
#include <string>

#include "common.hpp"

void set_options_type(rocksdb::Options &options, sm_idx_type type)
{
    if (type == K2I) {
        options.merge_operator.reset(new IDListOperator());
    }

    if (type == I2P) {
        options.merge_operator.reset(new PositionsMapOperator());
    }
}

void set_options_filter(rocksdb::Options &options)
{
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
}

void set_options_merge(rocksdb::Options &options)
{
    options.create_if_missing = true;
    options.IncreaseParallelism(4);
}

void open_filter(rocksdb::DB** db, const sm_config &conf, sm_idx_type type,
                 sm_idx_set set, int pid, bool ro)
{
    std::ostringstream rdb;
    rdb << conf.output_path << "/filter-" << sm::types[type] << "-"
        << sm::sets[set] << "." << pid << ".rdb";

    rocksdb::Status s;
    rocksdb::Options options;
    set_options_type(options, type);

    if (ro) {
        s = rocksdb::DB::OpenForReadOnly(options, rdb.str(), db);
    } else {
        set_options_filter(options);
        s = rocksdb::DB::Open(options, rdb.str(), db);
    }

    assert(s.ok());
}

void open_merge(rocksdb::DB** db, const sm_config &conf, sm_idx_type type,
                sm_idx_set set, bool ro)
{
    std::ostringstream rdb;
    rdb << conf.output_path << "/filter-" << sm::types[type] << "-"
        << sm::sets[set] << ".rdb";

    rocksdb::Status s;
    rocksdb::Options options;
    set_options_type(options, type);

    if (ro) {
        s = rocksdb::DB::OpenForReadOnly(options, rdb.str(), db);
    } else {
        set_options_merge(options);
        s = rocksdb::DB::Open(options, rdb.str(), db);
    }

    assert(s.ok());
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
