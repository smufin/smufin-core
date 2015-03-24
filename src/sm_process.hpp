#ifndef __SM_PROCESS_H__
#define __SM_PROCESS_H__

#include <sm_common.hpp>

#include <zlib.h>
#include <string>

void process_load(int pid, int lid);
void process_load_file(int pid, int lid, std::string file);
inline void process_load_sub(int pid, int lid, const char* sub, int len,
                             sm_read_kind kind, sm_bulk* bulks);

void process_incr(int sid, int num_loaders);
inline void process_incr_key(int sid, sm_key key, sm_read_kind kind);

#endif
