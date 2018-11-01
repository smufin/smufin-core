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

#include "index_format_plain.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

using std::cout;
using std::endl;
using std::string;

void index_format_plain::update(int fid, const sm_read *read, int pos,
                                char kmer[], sm_dir dir, sm_idx_set set)
{
    char buf[512] = {0};
    sprintf(buf, "%s %s", read->id, read->seq);

    _mutex[set].lock();
    auto result = _ids[set].insert(read->id);
    if (result.second == true) {
        _seq[set].insert(buf);
    }
    if (set == TM) {
        if (dir == DIR_A)
            _i2p[read->id].a[pos / 64] |= 1UL << (pos % 64);
        else
            _i2p[read->id].b[pos / 64] |= 1UL << (pos % 64);
    } else if (_k2i[set][kmer].size() <= _conf.max_filter_reads) {
        _k2i[set][kmer].insert(read->id);
    }
    _mutex[set].unlock();
}

bool index_format_plain::flush()
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

void index_format_plain::stats()
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

void index_format_plain::dump()
{
    std::vector<std::thread> threads;
    threads.push_back(std::thread(&index_format_plain::write_seq, this, NN));
    threads.push_back(std::thread(&index_format_plain::write_seq, this, TN));
    threads.push_back(std::thread(&index_format_plain::write_seq, this, TM));
    threads.push_back(std::thread(&index_format_plain::write_k2i, this, NN));
    threads.push_back(std::thread(&index_format_plain::write_k2i, this, TN));
    threads.push_back(std::thread(&index_format_plain::write_i2p, this, TM));
    cout << "Spawned " << threads.size() << " dump threads" << endl;
    for (auto& t: threads)
        t.join();
}

void index_format_plain::write_seq(sm_idx_set set)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path_filter << "/index-seq-" << sm::sets[set] << "."
         << _conf.pid << ".txt";
    ofs.open(file.str(), std::ofstream::app);
    for (auto const &s: _seq[set]) {
        ofs << s << "\n";
    }
    ofs.close();
}

void index_format_plain::write_k2i(sm_idx_set set)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path_filter << "/index-k2i-" << sm::sets[set] << "."
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

void index_format_plain::write_i2p(sm_idx_set set)
{
    std::ofstream ofs;
    std::ostringstream file;
    file << _conf.output_path_filter << "/index-i2p-" << sm::sets[set] << "."
         << _conf.pid << ".txt";
    ofs.open(file.str());
    for (auto const &kv: _i2p) {
        const sm_pos_bitmap *p = &kv.second;
        ofs << kv.first;
        for (int i = 0; i < POS_LEN; i++)
            ofs << " " << p->a[i];
        for (int i = 0; i < POS_LEN; i++)
            ofs << " " << p->b[i];
        ofs << "\n";
    }
    ofs.close();
}
