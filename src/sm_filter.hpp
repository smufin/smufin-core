#ifndef __SM_FILTER_H__
#define __SM_FILTER_H__

#include <zlib.h>
#include <string>
#include <kseq.h>

KSEQ_INIT(gzFile, gzread);

const char alpha[] = "ACGT";

void filter(int pid, int fid);
void filter_file(int pid, int fid, std::string file);

void filter_normal(int pid, int fid, kseq_t *seq, const char *sub, int len);
void filter_cancer(int pid, int fid, kseq_t *seq, const char *sub, int len);

void get_value(int pid, int fid, char kmer[], sm_tally *tally);

inline void filter_all(int pid, int fid, kseq_t *seq, int pos, bool rev,
                       char kmer[], sm_set set);
inline void filter_branch(int pid, int fid, kseq_t *seq, int pos, bool rev,
                          char kmer[], sm_set set);
inline void filter_kmer(kseq_t *seq, int pos, bool rev, char kmer[],
                        uint32_t nc,uint32_t tc, uint32_t nsum, uint32_t tsum,
                        sm_set set);

#endif
