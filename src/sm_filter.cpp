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
    std::size_t found = file.find("_N_");
    if (found != std::string::npos) {
        return;
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
                filter_sub(pid, fid, seq, &seq->seq.s[p], n);
                p += n;
            }
            p++;
        }

        n = 80 - p;
        if (n > 0) {
            filter_sub(pid, fid, seq, &seq->seq.s[p], n);
        }
    }

    kseq_destroy(seq);
    gzclose(in);
}

void filter_sub(int pid, int fid, kseq_t *seq, const char *sub, int len)
{
    char kmer[31];
    if (len >= KMER_LEN) {
        for (int i = 0; i <= len - KMER_LEN; i++) {
            strncpy(kmer, &sub[i], KMER_LEN);
            kmer[30] = '\0';

            uint32_t m = 0;
            memcpy(&m, kmer, MAP_LEN);
            hash_4c_map(m);

            if (map_l1[m] != pid)
                continue;
            int sid = map_l2[m];
            sm_key key = strtob4(kmer);

            sm_value value = tables[sid][key];
            uint32_t cnr = value.first;
            uint32_t ctr = value.second;
            char buf[256] = {0};
            if (ctr >= 4 && cnr == 0) {
                sprintf(buf, "@%s\n%s\n+\n%s",
                        seq->name.s, seq->seq.s, seq->qual.s);
                filter_mutex.lock();
                filter_reads.insert(buf);
                filter_mutex.unlock();
                break;
            }
        }
    }
}
