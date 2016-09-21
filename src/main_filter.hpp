#ifndef __SM_MAIN_FILTER_H__
#define __SM_MAIN_FILTER_H__

#include <fstream>

int main(int argc, char *argv[]);
void display_usage();

void reset_input_queue(std::ifstream &input_file);
void init();
void free_tables();
void free_table(int sid);
void rebuild_tables(int pid);
void rebuild_table(int pid, int sid);

void sm_process(int pid, int num_loaders, int num_storers);
void sm_filter(int pid, int num_filters);
void sm_stats(int num_storers);

void sm_write_fastq(int set, int pid);
void sm_write_k2i(int set, int pid);
void sm_write_i2p(int set, int pid);

#endif
