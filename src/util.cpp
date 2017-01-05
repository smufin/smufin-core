#include "util.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>

using std::cout;
using std::endl;
using std::string;

// Generic and timed thread spawn, to avoid duplicating the same boilerplate
// for different stages/steps.
template<typename T>
void spawn(string name, std::function<void(T)> func, std::vector<T> list)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    std::vector<std::thread> threads;
    for (auto& i: list)
        threads.push_back(std::thread(func, i));
    cout << "Spawned " << threads.size() << " " << name << " threads" << endl;
    for (auto& thread: threads)
        thread.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Spawn " << name << " time: " << time.count() << endl;
}

// Specialization to spawn N threads, passing a list of sequential integers:
// [0..N).
void spawn(string name, std::function<void(int)> func, int n)
{
    std::vector<int> list;
    for (int i = 0; i < n; i++)
        list.push_back(i);
    spawn<int>(name, func, list);
}

// Estimate size in GB of a sparsehash table of n elements, with keys of size
// k and values of size v.
float estimate_sparse(uint64_t n, size_t k, size_t v)
{
    return (n * (k + v + 1)) / 1024 / 1024 / 1024;
}
