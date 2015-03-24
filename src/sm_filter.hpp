#ifndef __SM_FILTER_H__
#define __SM_FILTER_H__

#include <zlib.h>
#include <string>
#include <kseq.h>

KSEQ_INIT(gzFile, gzread);

void filter(int pid, int fid);
void filter_file(int pid, int fid, std::string file);
void filter_sub(int pid, int fid, kseq_t *seq, const char *sub, int len);

#endif
