#ifndef __SM_STANDALONE_H__
#define __SM_STANDALONE_H__

#include <fstream>

int main(int argc, char *argv[]);
void display_usage();

void reset_input_queue(std::ifstream &input_file);
void sm_process(int pid, int num_loaders, int num_storers);
void sm_filter(int pid, int num_filters);
void sm_stats(int num_storers);

#endif
