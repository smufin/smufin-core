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

#include "db.hpp"

#include <iostream>
#include <sstream>
#include <string>

#include <rocksdb/cache.h>
#include <rocksdb/filter_policy.h>
#include <rocksdb/table.h>
#include <rocksdb/utilities/options_util.h>

#include "common.hpp"

using std::cout;
using std::endl;
using std::string;

void set_options_type(rocksdb::ColumnFamilyOptions &options, sm_idx_type type)
{
    if (type == K2I) {
        options.merge_operator.reset(new IDListOperator());
    }

    if (type == I2P) {
        options.merge_operator.reset(new PositionsMapOperator());
    }
}

void open_index_part_load(const sm_config &conf, sm_idx_type type,
                          sm_idx_set set, int pid, int iid,
                          rdb_handle &rdb)
{
    std::ostringstream conf_file, path;
    conf_file << conf.data_path << "/rocks/filter.conf";
    path << conf.output_path_filter << "/index-" << sm::types[type] << "-"
         << sm::sets[set] << "." << pid << "-" << iid << ".rdb";
    open_index(conf, type, path.str(), conf_file.str(), rdb);
}

void open_index_part_iter(const sm_config &conf, sm_idx_type type,
                          sm_idx_set set, int pid, int iid,
                          rdb_handle &rdb)
{
    std::ostringstream path;
    path << conf.output_path_filter << "/index-" << sm::types[type] << "-"
         << sm::sets[set] << "." << pid << "-" << iid << ".rdb";
    open_index_ro(conf, type, path.str(), rdb);
}

void open_index_full_load(const sm_config &conf, sm_idx_type type,
                          sm_idx_set set, rdb_handle &rdb)
{
    std::ostringstream conf_file, path;
    conf_file << conf.data_path << "/rocks/merge.conf";
    path << conf.output_path_merge << "/index-" << sm::types[type] << "-"
         << sm::sets[set] << ".rdb";
    open_index(conf, type, path.str(), conf_file.str(), rdb);
}

void open_index_full_iter(const sm_config &conf, sm_idx_type type,
                          sm_idx_set set, rdb_handle &rdb)
{
    std::ostringstream path;
    path << conf.output_path_merge << "/index-" << sm::types[type] << "-"
         << sm::sets[set] << ".rdb";
    open_index_ro(conf, type, path.str(), rdb);
}

// Create and open RocksDB index, pointed by `path', with read-write access.
// Loads RocksDB options using the configuration file pointed by `conf_file'.
void open_index(const sm_config &conf, sm_idx_type type,
                const std::string &path, const std::string &conf_file,
                rdb_handle &rdb)
{
    rocksdb::Status s;
    rocksdb::DBOptions db_options;
    std::vector<rocksdb::ColumnFamilyDescriptor> cf_descs;

    rocksdb::Env* env = rocksdb::Env::Default();
    env->SetBackgroundThreads(conf.num_threads_high, Env::Priority::HIGH);
    env->SetBackgroundThreads(conf.num_threads_low, Env::Priority::LOW);

    s = rocksdb::LoadOptionsFromFile(conf_file, env, &db_options, &cf_descs);
    if (!s.ok()) {
        cout << "Failed to load RocksDB options: " << conf_file << endl;
        exit(1);
    }

    set_options_type(cf_descs[0].options, type);
    db_options.error_if_exists = true;

    rocksdb::BlockBasedTableOptions t_options;
    t_options.block_cache = rocksdb::NewLRUCache(conf.block_cache_size);
    t_options.block_size = conf.block_size;
    t_options.pin_l0_filter_and_index_blocks_in_cache = true;
    auto t_factory = rocksdb::NewBlockBasedTableFactory(t_options);
    cf_descs[0].options.table_factory.reset(t_factory);

    s = rocksdb::DB::Open(db_options, path, cf_descs, &rdb.cfs, &rdb.db);
    if (!s.ok()) {
        cout << "Failed to open RocksDB database: " << path << endl;
        exit(1);
    }
}

// Open an existing RocksDB index, pointed by `path', with read-only access.
void open_index_ro(const sm_config &conf, sm_idx_type type,
                   const std::string &path, rdb_handle &rdb)
{
    rocksdb::Status s;
    rocksdb::DBOptions db_options;
    std::vector<rocksdb::ColumnFamilyDescriptor> cf_descs;
    rocksdb::Env* env = rocksdb::Env::Default();

    s = rocksdb::LoadLatestOptions(path, env, &db_options, &cf_descs);
    if (!s.ok()) {
        cout << "Failed to load latest RocksDB options: " << path << endl;
        exit(1);
    }

    set_options_type(cf_descs[0].options, type);
    db_options.create_if_missing = false;
    db_options.wal_dir = path;

    s = rocksdb::DB::OpenForReadOnly(db_options, path, cf_descs, &rdb.cfs,
                                     &rdb.db);
    if (!s.ok()) {
        cout << "Failed to open RocksDB database: " << path << endl;
        exit(1);
    }
}

// Open fully merged index, optimizing for random reading. Unlike
// open_index_full_iter, open_index_full_read can make use of bigger caches
// and bloom filters to speed up non-sequential retrieval.
void open_index_full_read(const sm_config &conf, sm_idx_type type,
                          sm_idx_set set, rdb_handle &rdb)
{
    std::ostringstream conf_file, path;
    conf_file << conf.data_path << "/rocks/group.conf";
    path << conf.output_path_merge << "/index-" << sm::types[type] << "-"
         << sm::sets[set] << ".rdb";

    rocksdb::Status s;
    rocksdb::DBOptions db_options;
    std::vector<rocksdb::ColumnFamilyDescriptor> cf_descs;
    rocksdb::Env* env = rocksdb::Env::Default();

    s = rocksdb::LoadOptionsFromFile(conf_file.str(), env, &db_options,
                                     &cf_descs);
    if (!s.ok()) {
        cout << "Failed to load RocksDB options: " << conf_file.str() << endl;
        exit(1);
    }

    set_options_type(cf_descs[0].options, type);
    db_options.create_if_missing = false;
    cf_descs[0].options.disable_auto_compactions = true;

    rocksdb::BlockBasedTableOptions t_options;
    t_options.block_cache = rocksdb::NewLRUCache(conf.block_cache_size);
    t_options.block_size = conf.block_size;
    t_options.pin_l0_filter_and_index_blocks_in_cache = true;
    t_options.filter_policy.reset(rocksdb::NewBloomFilterPolicy(8));
    auto t_factory = rocksdb::NewBlockBasedTableFactory(t_options);
    cf_descs[0].options.table_factory.reset(t_factory);

    s = rocksdb::DB::OpenForReadOnly(db_options, path.str(), cf_descs,
                                     &rdb.cfs, &rdb.db);
    if (!s.ok()) {
        cout << "Failed to open RocksDB database: " << path.str() << endl;
        exit(1);
    }
}

void encode_pos(const sm_pos_bitmap &p, std::string &s)
{
    std::stringstream e;
    for (int i = 0; i < POS_LEN; i++)
        e << std::hex << p.a[i] << " ";
    for (int i = 0; i < POS_LEN; i++)
        e << std::hex << p.b[i] << " ";
    s = e.str();
}

sm_pos_bitmap decode_pos(const std::string &s)
{
    sm_pos_bitmap p;
    std::istringstream in(s);
    for (int i = 0; i < POS_LEN; i++)
        in >> std::hex >> p.a[i];
    for (int i = 0; i < POS_LEN; i++)
        in >> std::hex >> p.b[i];
    return p;
}
