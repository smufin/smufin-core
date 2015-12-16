#include <sm_common.hpp>
#include <sm_filter.hpp>

#include <string>
#include <iostream>
#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::string;

void filter(int pid, int fid)
{
    string file;
    while (input_count > 0) {
        while (input_queue.try_dequeue(file)) {
            filter_file(pid, fid, file);
            input_count--;
        }
    }
}

void filter_file(int pid, int fid, string file)
{
    // Identify read kind from file name.
    sm_read_kind kind = NORMAL_READ;
    std::size_t found = file.find("_T_");
    if (found != std::string::npos)
        kind = CANCER_READ;

    int len;
    gzFile in = gzopen(file.c_str(), "rb");

    kseq_t *seq = kseq_init(in);
    while ((len = kseq_read(seq)) >= 0) {
        if (lq_count(seq->qual.s) > 8)
            continue;

        int p = 0;
        int n = 80;
        char *ps;

        while ((ps = (char*) memchr(&seq->seq.s[p], 'N', 80 - p)) != NULL) {
            n = ps - &seq->seq.s[p];
            if (n > 0) {
                if (kind == CANCER_READ)
                    filter_cancer(pid, fid, seq, &seq->seq.s[p], n);
                else
                    filter_normal(pid, fid, seq, &seq->seq.s[p], n);
                p += n;
            }
            p++;
        }

        n = 80 - p;
        if (n > 0) {
            if (kind == CANCER_READ)
                filter_cancer(pid, fid, seq, &seq->seq.s[p], n);
            else
                filter_normal(pid, fid, seq, &seq->seq.s[p], n);
        }
    }

    kseq_destroy(seq);
    gzclose(in);
}

void filter_normal(int pid, int fid, kseq_t *seq, const char *sub, int len)
{
    if (len < KMER_LEN)
        return;

    char kmer[KMER_LEN + 1];
    for (int i = 0; i <= len - KMER_LEN; i++) {
        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';
        filter_all(pid, fid, seq, i, false, kmer, NN);

        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';
        krevcomp(kmer);
        filter_all(pid, fid, seq, i, true, kmer, NN);
    }
}

void filter_cancer(int pid, int fid, kseq_t *seq, const char *sub, int len)
{
    if (len < KMER_LEN)
        return;

    char kmer[KMER_LEN + 1];
    for (int i = 0; i <= len - KMER_LEN; i++) {
        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';
        filter_branch(pid, fid, seq, i, false, kmer, TM);
        filter_all(pid, fid, seq, i, false, kmer, TN);

        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';
        krevcomp(kmer);
        filter_branch(pid, fid, seq, i, true, kmer, TM);
        filter_all(pid, fid, seq, i, true, kmer, TN);
    }
}

int get_value(int pid, int fid, char kmer[], sm_table::const_iterator *it)
{
    char last = kmer[KMER_LEN - 1];
    kmer[KMER_LEN - 1] = '\0';

    uint32_t m = 0;
    memcpy(&m, &kmer[1], MAP_LEN);
    hash_4c_map(m);

    if (map_l1[m] != pid)
        return -1;
    int sid = map_l2[m];
    sm_key key = strtob4(&kmer[1]);

    *it = tables[sid].find(key);
    if (*it == tables[sid].end())
        return -1;

    kmer[KMER_LEN - 1] = last;
    return 0;
}

void filter_branch(int pid, int fid, kseq_t *seq, int pos, bool rev,
                   char kmer[], sm_set set)
{
    sm_table::const_iterator it;
    if (get_value(pid, fid, kmer, &it) != 0)
        return;
    int f = code[kmer[0]] - '0';
    int l = code[kmer[KMER_LEN - 1]] - '0';
    uint32_t nc = it->second.v[f][l][NORMAL_READ];
    uint32_t tc = it->second.v[f][l][CANCER_READ];
    uint32_t nsum = 0;
    uint32_t tsum = 0;
    for (l = 0; l < 4; l++) {
        nsum += it->second.v[f][l][NORMAL_READ];
        tsum += it->second.v[f][l][CANCER_READ];
    }
    filter_kmer(seq, pos, rev, kmer, nc, tc, nsum, tsum, set);
}

void filter_all(int pid, int fid, kseq_t *seq, int pos, bool rev,
                char kmer[], sm_set set)
{
    sm_table::const_iterator it;
    if (get_value(pid, fid, kmer, &it) != 0)
        return;
    for (int f = 0; f < 4; f++) {
        kmer[0] = alpha[f];
        uint32_t nsum = 0;
        uint32_t tsum = 0;
        for (int l = 0; l < 4; l++) {
            nsum += it->second.v[f][l][NORMAL_READ];
            tsum += it->second.v[f][l][CANCER_READ];
        }
        for (int l = 0; l < 4; l++) {
            kmer[KMER_LEN - 1] = alpha[l];
            uint32_t nc = it->second.v[f][l][NORMAL_READ];
            uint32_t tc = it->second.v[f][l][CANCER_READ];
            filter_kmer(seq, pos, rev, kmer, nc, tc, nsum, tsum, set);
        }
    }
}

void filter_kmer(kseq_t *seq, int pos, bool rev, char kmer[], uint32_t nc,
                 uint32_t tc, uint32_t nsum, uint32_t tsum, sm_set set)
{
    if (tc >= 4 && nc == 0) {
        char buf[256] = {0};
        sprintf(buf, "@%s\n%s\n+\n%s", seq->name.s, seq->seq.s, seq->qual.s);
        filter_mutex[set].lock();
        filter_reads[set].insert(buf);
        if (!rev)
            filter_i2p[set][seq->name.s].first.push_back((uint8_t) pos);
        else
            filter_i2p[set][seq->name.s].second.push_back((uint8_t) pos);
        filter_k2i[set][kmer].insert(seq->name.s);
        filter_mutex[set].unlock();
    }
}
