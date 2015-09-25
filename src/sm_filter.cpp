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
    if (found != std::string::npos) {
        kind = CANCER_READ;
    }

    int len;
    gzFile in = gzopen(file.c_str(), "rb");

    kseq_t *seq = kseq_init(in);
    while ((len = kseq_read(seq)) >= 0) {
        if (lq_count(seq->qual.s) > 8) {
            continue;
        }

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

    const char vars[] = "ACGT";

    char kmer[KMER_LEN + 1];
    for (int i = 0; i <= len - KMER_LEN; i++) {
        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';
        for (int j = 0; j < 4; j++) {
            kmer[0] = vars[j];

            uint32_t normal_parent = 0;
            uint32_t tumour_parent = 0;
            uint32_t normal_counts[4] = {0};
            uint32_t tumour_counts[4] = {0};

            for (int k = 0; k < 4; k++) {
                kmer[KMER_LEN - 1] = vars[k];

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

                uint32_t cnr = it->second.first;
                uint32_t ctr = it->second.second;

                normal_parent += cnr;
                tumour_parent += ctr;
                normal_counts[k] = cnr;
                tumour_counts[k] = ctr;
            }

            for (int k = 0; k < 4; k++) {
                kmer[KMER_LEN - 1] = vars[k];

                uint32_t cnr = normal_counts[k];
                uint32_t ctr = tumour_counts[k];

                if (ctr >= 4 && cnr == 0) {
                    char buf[256] = {0};
                    sprintf(buf, "@%s\n%s\n+\n%s",
                            seq->name.s, seq->seq.s, seq->qual.s);
                    filter_nn_mutex.lock();
                    filter_nn_reads.insert(buf);
                    filter_nn_mutex.unlock();
                }
            }
        }
    }
}

void filter_cancer(int pid, int fid, kseq_t *seq, const char *sub, int len)
{
    if (len < KMER_LEN)
        return;

    const char vars[] = "ACGT";

    char kmer[KMER_LEN + 1];

    for (int i = 0; i <= len - KMER_LEN; i++) {
        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';

        char last = kmer[KMER_LEN - 1];
        uint32_t normal_parent = 0;
        uint32_t tumour_parent = 0;
        uint32_t normal_counts[4] = {0};
        uint32_t tumour_counts[4] = {0};

        for (int k = 0; k < 4; k++) {
            kmer[KMER_LEN - 1] = vars[k];

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

            uint32_t cnr = it->second.first;
            uint32_t ctr = it->second.second;

            normal_parent += cnr;
            tumour_parent += ctr;
            normal_counts[k] = cnr;
            tumour_counts[k] = ctr;
        }

        kmer[KMER_LEN - 1] = last;
        uint32_t cnr = normal_counts[code[last] - '0'];
        uint32_t ctr = tumour_counts[code[last] - '0'];

        if (ctr >= 4 && cnr == 0) {
            char buf[256] = {0};
            sprintf(buf, "@%s\n%s\n+\n%s",
                    seq->name.s, seq->seq.s, seq->qual.s);
            filter_tm_mutex.lock();
            filter_tm_reads.insert(buf);
            filter_tm_mutex.unlock();
        }
    }

    for (int i = 0; i <= len - KMER_LEN; i++) {
        strncpy(kmer, &sub[i], KMER_LEN);
        kmer[KMER_LEN] = '\0';
        for (int j = 0; j < 4; j++) {
            kmer[0] = vars[j];

            uint32_t normal_parent = 0;
            uint32_t tumour_parent = 0;
            uint32_t normal_counts[4] = {0};
            uint32_t tumour_counts[4] = {0};

            for (int k = 0; k < 4; k++) {
                kmer[KMER_LEN - 1] = vars[k];

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

                uint32_t cnr = it->second.first;
                uint32_t ctr = it->second.second;

                normal_counts[k] = cnr;
                tumour_counts[k] = ctr;
                normal_parent += cnr;
                tumour_parent += ctr;
            }

            for (int k = 0; k < 4; k++) {
                kmer[KMER_LEN - 1] = vars[k];

                uint32_t cnr = normal_counts[k];
                uint32_t ctr = tumour_counts[k];

                if (ctr >= 4 && cnr == 0) {
                    char buf[256] = {0};
                    sprintf(buf, "@%s\n%s\n+\n%s",
                            seq->name.s, seq->seq.s, seq->qual.s);
                    filter_tn_mutex.lock();
                    filter_tn_reads.insert(buf);
                    filter_tn_mutex.unlock();
                }
            }
        }
    }
}
