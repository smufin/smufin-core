#ifndef __SM_MAIN_FILTER_H__
#define __SM_MAIN_FILTER_H__

#include <fstream>

int main(int argc, char *argv[]);
void display_usage();

void reset_input_queue(std::ifstream &input_file);
void free_caches();

void sm_process(int pid, int num_loaders, int num_storers);
void sm_filter(int pid, int num_filters);
void sm_stats(int num_storers);

void sm_write_fastq(int set);
void sm_write_k2i(int set);
void sm_write_i2p(int set);

#endif
