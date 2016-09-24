#ifndef __SM_MAIN_H__
#define __SM_MAIN_H__

#include <common.hpp>
#include <stage.hpp>

int main(int argc, char *argv[]);
void display_usage();

template<typename T> stage* create_stage(const sm_config &conf);

#endif
