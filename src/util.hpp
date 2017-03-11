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

#endif
