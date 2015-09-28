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
        filter_tree(pid, fid, seq, kmer, NN_WAY);
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

        char last = kmer[KMER_LEN - 1];
        uint32_t nsum = 0;
        uint32_t tsum = 0;
        uint32_t narr[4] = {0};
        uint32_t tarr[4] = {0};
        get_branch(pid, fid, kmer, narr, tarr, &nsum, &tsum);

        kmer[KMER_LEN - 1] = last;
        uint32_t nc = narr[code[last] - '0'];
        uint32_t tc = tarr[code[last] - '0'];
        filter_kmer(seq, kmer, nc, tc, nsum, tsum, TM_WAY);

        filter_tree(pid, fid, seq, kmer, TN_WAY);
    }
}

// Doesn't guarantee integrity of kmer's last character.
void get_branch(int pid, int fid, char kmer[], uint32_t narr[],
                uint32_t tarr[], uint32_t *nsum, uint32_t *tsum)
{
    for (int i = 0; i < 4; i++) {
        kmer[KMER_LEN - 1] = alpha[i];

        uint32_t m = 0;
        memcpy(&m, &kmer[0], MAP_LEN);
        hash_4c_map(m);

        if (map_l1[m] != pid)
            continue;
        int sid = map_l2[m];
        sm_key key = strtob4(kmer);

        sm_table::const_iterator it = tables[sid].find(key);
        if (it == tables[sid].end())
            continue;

        uint32_t nc = it->second.first;
        uint32_t tc = it->second.second;

        *nsum += nc;
        *tsum += tc;
        narr[i] = nc;
        tarr[i] = tc;
    }
}

void filter_tree(int pid, int fid, kseq_t *seq, char kmer[], sm_way way)
{
    for (int i = 0; i < 4; i++) {
        kmer[0] = alpha[i];
        uint32_t nsum = 0;
        uint32_t tsum = 0;
        uint32_t narr[4] = {0};
        uint32_t tarr[4] = {0};
        get_branch(pid, fid, kmer, narr, tarr, &nsum, &tsum);
        for (int j = 0; j < 4; j++) {
            kmer[KMER_LEN - 1] = alpha[j];
            filter_kmer(seq, kmer, narr[j], tarr[j], nsum, tsum, way);
        }
    }
}

void filter_kmer(kseq_t *seq, char kmer[], uint32_t nc, uint32_t tc,
                 uint32_t nsum, uint32_t tsum, sm_way way)
{
    if (tc >= 4 && nc == 0) {
        char buf[256] = {0};
        sprintf(buf, "@%s\n%s\n+\n%s", seq->name.s, seq->seq.s, seq->qual.s);
        filter_mutex[way].lock();
        filter_reads[way].insert(buf);
        filter_mutex[way].unlock();
    }
}
