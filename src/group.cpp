#include "group.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::string;

void rrevcomp(char read[], int len)
{
    const char comp[] = "-------------------------------------------"
                        "----------------------T-G---C------N-----A";
    int c, i, j;
    for (i = 0, j = len - 1; i < j; i++, j--) {
        c = read[i];
        read[i] = comp[read[j]];
        read[j] = comp[c];
    }
}

void get_positions_a(uint64_t bitmap[2], std::vector<int> *pos)
{
    for (int i = 0; i < 2; i++) {
        unsigned long tmp = bitmap[i];
        int offset = i * 64;
        while (tmp > 0) {
            int p = __builtin_ffsl(tmp) - 1;
            tmp &= (tmp - 1);
            pos->push_back(p + offset);
        }
    }
}

void get_positions_b(uint64_t bitmap[2], std::vector<int> *pos, int len)
{
    for (int i = 0; i < 2; i++) {
        unsigned long tmp = bitmap[i];
        int offset = i * 64;
        while (tmp > 0) {
            int p = __builtin_ffsl(tmp) - 1;
            tmp &= (tmp - 1);
            pos->push_back(len - KMER_LEN - (p + offset));
        }
    }
}

bool match_window(std::vector<int> pos)
{
    if (pos.size() < WMIN) {
        return false;
    }

    for (int i = 0; i <= pos.size() - WMIN; i++) {
        std::vector<int> sub(pos.begin() + i, pos.begin() + i + WMIN);
        if (sub.back() - sub.front() < WLEN) {
            return true;
        }
    }

    return false;
}

void select_candidate(int gid, string sid, string seq, std::vector<int>& pos,
                      int dir)
{
    (*l2p[gid])[sid][dir] = pos;
    // TODO: Switch to unordered_set instead of vector?
    std::vector<string> kmers;
    for (int p: pos) {
        kmers.push_back(seq.substr(p, KMER_LEN));
    }
    (*l2k[gid])[sid][dir] = kmers;
}

void populate(int pid, int gid)
{
    const char comp_code[] = "ab";
    const char kind_code[] = "nt";

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    std::ofstream ofs;
    ofs.open("groups." + std::to_string(pid) + "-" + std::to_string(gid) + ".json");

    int num_leads = 0;
    ofs << "{";
    bool first_group = true;
    for (l2k_table::const_iterator it = l2k[gid]->begin(); it != l2k[gid]->end(); ++it) {
        string lid = it->first;
        index_count keep;
        index_count drop;

        for (int i = 0; i < 2; i++) {
            for (string kmer: it->second[i]) {
                keep[i][kmer] = 0;
                drop[i][kmer] = 0;
            }
        }

        populate_index(gid, lid, it->second[0], NN, keep, drop);
        populate_index(gid, lid, it->second[0], TN, keep, drop);
        populate_index(gid, lid, it->second[1], NN, keep, drop);
        populate_index(gid, lid, it->second[1], TN, keep, drop);

        string read;
        rocksdb::Status status;
        status = i2r[TM]->Get(rocksdb::ReadOptions(), lid, &read);
        if (!status.ok())
            continue;

        if (!first_group)
            ofs << ",";
        first_group = false;

        ofs << "\"" << lid << "\":{";
        ofs  << "\"lead\":["
             << "\"" << lid << "\","
             << "\"" << read << "\""
             << "],";

        for (int i = 0; i < 2; i++) {
            ofs << "\"pos-" << comp_code[i] << "\":[";
            bool first_pos = true;
            for (int p: (*l2p[gid])[lid][i]) {
                if (!first_pos)
                    ofs << ",";
                first_pos = false;
                ofs << p;
            }
            ofs << "],";

            ofs << "\"kmers-" << comp_code[i] << "\":[";
            bool first_kmer = true;
            for (string kmer: it->second[i]) {
                int kept_n = keep[0][kmer];
                int dropped_n = drop[0][kmer];
                int kept_t = keep[1][kmer];
                int dropped_t = drop[1][kmer];
                if (!first_kmer)
                    ofs << ",";
                first_kmer = false;
                ofs << "[\"" << kmer << "\"," << kept_n << "," << kept_t << ","
                     << dropped_n << "," << dropped_t << "]";
            }
            ofs << "],";
        }

        for (int i = 0; i < 2; i++) {
            bool first_read = true;
            ofs << "\"reads-" << kind_code[i] << "\":[";
            for (string sid: (*l2i[gid])[lid][i]) {
                if (!first_read)
                    ofs << ",";
                first_read = false;
                read = "";
                status = i2r[i]->Get(rocksdb::ReadOptions(), lid, &read);
                if (!status.ok())
                    continue;
                ofs << "[\"" << sid << "\",\"" << read << "\"]";
            }
            ofs << ( i == 1 ? "]" : "]," );
        }

        ofs << "}";

        if (num_leads % 100 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "P: " << gid << " " << num_leads << " "
                 << time.count() << "\n";
            start = std::chrono::system_clock::now();
        }
        num_leads++;
    }
    ofs << "}";
}

void populate_index(int gid, string& lid, const std::vector<string>& kmers,
                    int kind, index_count& keep, index_count& drop)
{
    for (string kmer: kmers) {
        string list;
        rocksdb::Status status;
        status = k2i[kind]->Get(rocksdb::ReadOptions(), kmer, &list);
        if (!status.ok())
            continue;

        std::unordered_set<string> sids;
        boost::split(sids, list, boost::is_any_of(" "));

        if (sids.size() > DROP) {
            drop[kind][kmer] += sids.size();
            continue;
        }

        (*l2i[gid])[lid][kind].insert(sids.begin(), sids.end());
        keep[kind][kmer] += sids.size();
    }
}
