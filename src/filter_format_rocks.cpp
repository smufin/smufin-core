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
    rocksdb::Options opts;
    rocksdb::Status status;

    for (int i = 0; i < NUM_SETS; i++) {
        opts = get_rocks_options("seq");
        std::ostringstream path;
        path << _conf.output_path << "/filter-seq" << "-" << sm::sets[i]
             << "." << _conf.pid << ".rdb";
        status = rocksdb::DB::Open(opts, path.str(), &_seq[i]);
        assert(status.ok());
    }

    for (int i = 0; i < NUM_SETS; i++) {
        opts = get_rocks_options("k2i");
        std::ostringstream path;
        path << _conf.output_path << "/filter-k2i" << "-" << sm::sets[i]
             << "." << _conf.pid << ".rdb";
        status = rocksdb::DB::Open(opts, path.str(), &_k2i[i]);
        assert(status.ok());
    }

    opts = get_rocks_options("i2p");
    std::ostringstream path;
    path << _conf.output_path << "/filter-i2p-tm." << _conf.pid << ".rdb";
    status = rocksdb::DB::Open(opts, path.str(), &_i2p);
    assert(status.ok());
}

void filter_format_rocks::update(kseq_t *seq, int pos, bool rev, char kmer[],
                                 sm_set set)
{
    string sid = seq->name.s;
    rocksdb::WriteOptions w_options;
    w_options.disableWAL = true;

    _seq[set]->Put(w_options, sid, seq->seq.s);

    // TODO: Honour _conf.max_k2i_reads
    std::stringstream ss;
    ss << sid << " ";
    _k2i[set]->Merge(w_options, kmer, ss.str());

    if (set == TM) {
        sm_pos_bitmap p;
        string serialized;

        if (!rev)
            p.a[pos / 64] |= 1UL << (pos % 64);
        else
            p.b[pos / 64] |= 1UL << (pos % 64);

        encode_pos(serialized, p);
        _i2p->Merge(w_options, sid, serialized);
    }
}

void filter_format_rocks::dump()
{
    std::vector<rocksdb::DB*> list;
    for (int i = 0; i < NUM_SETS; i++) {
        list.push_back(_seq[i]);
        list.push_back(_k2i[i]);
    }
    list.push_back(_i2p);
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

    uint64_t nn, tn, tm;
    rocksdb::Iterator* it;

    nn = tn = tm = 0;
    it = _seq[NN]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) nn++;
    it = _seq[TN]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tn++;
    it = _seq[TM]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tm++;
    cout << "Size SEQ: " << nn << " " << tn << " " << tm << endl;

    nn = tn = tm = 0;
    it = _k2i[NN]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) nn++;
    it = _k2i[TN]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tn++;
    it = _k2i[TM]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tm++;
    cout << "Size K2I: " << nn << " " << tn << " " << tm << endl;

    tm = 0;
    it = _i2p->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tm++;
    cout << "Size I2P: " << tm << endl;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Time filter/stats: " << time.count() << endl;
}
