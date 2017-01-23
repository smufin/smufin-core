#include "filter_format_plain.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

using std::cout;
using std::endl;
using std::string;

void filter_format_plain::update(kseq_t *seq, int pos, bool rev, char kmer[],
                                 sm_idx_set set)
{
    char buf[512] = {0};
    sprintf(buf, "%s %s", seq->name.s, seq->seq.s);

    _mutex[set].lock();
    auto result = _ids[set].insert(seq->name.s);
    if (result.second == true) {
        _seq[set].insert(buf);
    }
    if (set == TM) {
        if (!rev)
            _i2p[seq->name.s].a[pos / 64] |= 1UL << (pos % 64);
        else
            _i2p[seq->name.s].b[pos / 64] |= 1UL << (pos % 64);
    } else if (_k2i[set][kmer].size() <= _conf.max_filter_reads) {
        _k2i[set][kmer].insert(seq->name.s);
    }
    _mutex[set].unlock();
}

bool filter_format_plain::flush()
{
    bool flushed = false;
    for (auto set: {NN, TN, TM}) {
        _mutex[set].lock();
        if (_seq[set].size() > 1000000) {
            write_seq(set);
            _seq[set].clear();
            _seq[set] = std::unordered_set<std::string>();
            flushed = true;
        }
        _mutex[set].unlock();
    }
    return flushed;
}

void filter_format_plain::stats()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    cout << "Size SEQ: " << _ids[NN].size() << " " << _ids[TN].size() << " "
         << _ids[TM].size() << endl;
    cout << "Size K2I: " << _k2i[NN].size() << " " << _k2i[TN].size() << endl;
    cout << "Size I2P: " << _i2p.size() << endl;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Time filter/stats: " << time.count() << endl;
}

void filter_format_plain::dump()
{
    std::vector<std::thread> threads;
    threads.push_back(std::thread(&filter_format_plain::write_seq, this, NN));
    threads.push_back(std::thread(&filter_format_plain::write_seq, this, TN));
    threads.push_back(std::thread(&filter_format_plain::write_seq, this, TM));
    threads.push_back(std::thread(&filter_format_plain::write_k2i, this, NN));
    threads.push_back(std::thread(&filter_format_plain::write_k2i, this, TN));
    threads.push_back(std::thread(&filter_format_plain::write_i2p, this, TM));
    cout << "Spawned " << threads.size() << " dump threads" << endl;
    for (auto& t: threads)
        t.join();
}

void filter_format_plain::write_seq(sm_idx_set set)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path << "/filter-seq-" << sm::sets[set] << "."
         << _conf.pid << ".txt";
    ofs.open(file.str(), std::ofstream::app);
    for (auto const &s: _seq[set]) {
        ofs << s << "\n";
    }
    ofs.close();
}

void filter_format_plain::write_k2i(sm_idx_set set)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path << "/filter-k2i-" << sm::sets[set] << "."
         << _conf.pid << ".txt";
    ofs.open(file.str());
    for (auto const &kv: _k2i[set]) {
        if (kv.second.size() > _conf.max_filter_reads)
            continue;
        ofs << kv.first << " " << kv.second.size();
        for (auto const &sid: kv.second) {
            ofs << " " << sid;
        }
        ofs << "\n";
    }
    ofs.close();
}

void filter_format_plain::write_i2p(sm_idx_set set)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path << "/filter-i2p-" << sm::sets[set] << "."
         << _conf.pid << ".txt";
    ofs.open(file.str());
    for (auto const &kv: _i2p) {
        const sm_pos_bitmap *p = &kv.second;
        ofs << kv.first << " ";
        ofs << p->a[0] << " " << p->a[1] << " ";
        ofs << p->b[0] << " " << p->b[1] << "\n";
    }
    ofs.close();
}
