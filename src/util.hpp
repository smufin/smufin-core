/*
 * Copyright © 2015-2019 Barcelona Supercomputing Center (BSC)
 *
 * This file is part of SMUFIN Core. SMUFIN Core is released under the SMUFIN
 * Public License, and may not be used except in compliance with it. This file
 * is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; see the SMUFIN Public License for more details. You should have
 * received a copy of the SMUFIN Public License along with this file. If not,
 * see <https://github.com/smufin/smufin-core/blob/master/COPYING>.
 *
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2019
 */

#ifndef __SM_UTIL_H__
#define __SM_UTIL_H__

#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

#include "common.hpp"

template<typename T>
void spawn(std::string name, std::function<void(T)> func, std::vector<T> list)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    std::vector<std::thread> threads;
    for (auto& i: list)
        threads.push_back(std::thread(func, i));
    std::cout << "Spawned " << threads.size() << " " << name << " threads"
              << std::endl;
    for (auto& thread: threads)
        thread.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    std::cout << "Time " << name << " spawn: " << time.count() << std::endl;
}

void spawn(std::string name, std::function<void(int)> func, int n);

void init_mapping(const sm_config &conf, int n1, int n2, int l1[], int l2[]);

float estimate_sparse(uint64_t n, size_t k, size_t v);

std::vector<std::string> expand_path(std::string path);

// Functions to read data from serialized sparsehash tables.
bool read_be32(FILE* fp, uint32_t* value);
bool read_be64(FILE* fp, uint64_t* value);
bool read_be(FILE* fp, uint64_t* value);

#endif
