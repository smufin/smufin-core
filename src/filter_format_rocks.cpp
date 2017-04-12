#include "filter_format_rocks.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

#include "db.hpp"
#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

filter_format_rocks::filter_format_rocks(const sm_config &conf)
    : filter_format(conf)
{
    rocksdb::Status status;
    rocksdb::Options opts;
    set_options_filter(opts);

    for (int fid = 0; fid < _conf.num_filters; fid++) {
        for (auto set: {NN, TN, TM}) {
            std::ostringstream path;
            path << _conf.output_path << "/filter-seq" << "-" << sm::sets[set]
                 << "." << _conf.pid << "-" << fid << ".rdb";
            status = rocksdb::DB::Open(opts, path.str(), &_seq[set][fid]);
            assert(status.ok());
        }

        for (auto set: {NN, TN}) {
            set_options_type(opts, K2I);
            std::ostringstream path;
            path << _conf.output_path << "/filter-k2i" << "-" << sm::sets[set]
                 << "." << _conf.pid << "-" << fid << ".rdb";
            status = rocksdb::DB::Open(opts, path.str(), &_k2i[set][fid]);
            assert(status.ok());
        }

        set_options_type(opts, I2P);
        std::ostringstream path;
        path << _conf.output_path << "/filter-i2p-tm." << _conf.pid << "-"
             << fid << ".rdb";
        status = rocksdb::DB::Open(opts, path.str(), &_i2p[fid]);
        assert(status.ok());
    }
}

void filter_format_rocks::update(int fid, const sm_read *read, int pos,
                                 bool rev, char kmer[], sm_idx_set set)
{
    string sid = read->id;
    rocksdb::WriteOptions w_options;
    w_options.disableWAL = true;

    _seq[set][fid]->Put(w_options, sid, read->seq);

    if (set == TM) {
        sm_pos_bitmap p;
        string serialized;

        if (!rev)
            p.a[pos / 64] |= 1UL << (pos % 64);
        else
            p.b[pos / 64] |= 1UL << (pos % 64);

        encode_pos(p, serialized);
        _i2p[fid]->Merge(w_options, sid, serialized);
    } else {
        // TODO: Honour _conf.max_filter_reads
        std::stringstream ss;
        ss << sid << " ";
        _k2i[set][fid]->Merge(w_options, kmer, ss.str());
    }
}

void filter_format_rocks::dump()
{
    std::vector<rocksdb::DB*> list;
    for (int fid = 0; fid < _conf.num_filters; fid++) {
        for (auto set: {NN, TN, TM})
            list.push_back(_seq[set][fid]);
        for (auto set: {NN, TN})
            list.push_back(_k2i[set][fid]);
        list.push_back(_i2p[fid]);
    }
    spawn<rocksdb::DB*>("compact", std::bind(&filter_format_rocks::compact,
                        this, std::placeholders::_1), list);
}

void filter_format_rocks::compact(rocksdb::DB* db)
{
    rocksdb::CompactRangeOptions c_options = rocksdb::CompactRangeOptions();
    db->CompactRange(c_options, nullptr, nullptr);
}

void filter_format_rocks::stats()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    rocksdb::Iterator* it;
    uint64_t seq_nn, seq_tn, seq_tm;
    uint64_t k2i_nn, k2i_tn;
    uint64_t i2p_tm;

    seq_nn = seq_tn = seq_tm = k2i_nn = k2i_tn = i2p_tm = 0;

    for (int fid = 0; fid < _conf.num_filters; fid++) {
        it = _seq[NN][fid]->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) seq_nn++;
        it = _seq[TN][fid]->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) seq_tn++;
        it = _seq[TM][fid]->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) seq_tm++;

        it = _k2i[NN][fid]->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) k2i_nn++;
        it = _k2i[TN][fid]->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) k2i_tn++;

        it = _i2p[fid]->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) i2p_tm++;
    }

    cout << "Size SEQ: " << seq_nn << " " << seq_tn << " " << seq_tm << endl;
    cout << "Size K2I: " << k2i_nn << " " << k2i_tn << endl;
    cout << "Size I2P: " << i2p_tm << endl;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Time filter/stats: " << time.count() << endl;
}
