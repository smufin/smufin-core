#include "filter_format_rocks.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

#include "db.hpp"

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

    _seq[set]->Put(rocksdb::WriteOptions(), sid, seq->seq.s);

    // TODO: Honour _conf.max_k2i_reads
    std::stringstream s;
    s << sid << " ";
    _k2i[set]->Merge(rocksdb::WriteOptions(), kmer, s.str());

    if (set == TM) {
        sm_pos_bitmap p;
        string serialized;

        if (!rev)
            p.a[pos / 64] |= 1UL << (pos % 64);
        else
            p.b[pos / 64] |= 1UL << (pos % 64);

        encode_pos(serialized, p);
        _i2p->Merge(rocksdb::WriteOptions(), sid, serialized);
    }
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
    cout << "Count SEQ: " << nn << " " << tn << " " << tm << endl;

    nn = tn = tm = 0;
    it = _k2i[NN]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) nn++;
    it = _k2i[TN]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tn++;
    it = _k2i[TM]->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tm++;
    cout << "Count K2I: " << nn << " " << tn << " " << tm << endl;

    tm = 0;
    it = _i2p->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) tm++;
    cout << "Count I2P: " << tm << endl;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Time filter/stats: " << time.count() << endl;
}
