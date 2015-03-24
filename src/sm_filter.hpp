#ifndef __SM_FILTER_H__
#define __SM_FILTER_H__

#include <zlib.h>
#include <string>

void filter(int pid, int fid);
void filter_file(int pid, int fid, std::string file);
void filter_sub(int pid, int fid, const char* sub, int len);

#endif
