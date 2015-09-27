#ifndef __SM_FILTER_H__
#define __SM_FILTER_H__

#include <zlib.h>
#include <string>
#include <kseq.h>

KSEQ_INIT(gzFile, gzread);

const char alpha[] = "ACGT";

enum sm_way {
    NN_WAY, TN_WAY, TM_WAY
};

void filter(int pid, int fid);
void filter_file(int pid, int fid, std::string file);

void filter_normal(int pid, int fid, kseq_t *seq, const char *sub, int len);
void filter_cancer(int pid, int fid, kseq_t *seq, const char *sub, int len);

inline void get_branch(int pid, int fid, char kmer[], uint32_t narr[],
                       uint32_t tarr[], uint32_t *nsum, uint32_t *tsum);

inline void filter_tree(int pid, int fid, kseq_t *seq, char kmer[], sm_way way);
inline void filter_kmer(kseq_t *seq, char kmer[], uint32_t nc, uint32_t tc,
                        uint32_t nsum, uint32_t tsum, sm_way way);

#endif
