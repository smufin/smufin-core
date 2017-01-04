#include "common.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>

using std::cout;
using std::endl;
using std::string;

// String to integer conversion. Parses a null-terminated string, interpreting
// its content as an integral number in base 4.
uint64_t strtob4(const char *str)
{
    register uint64_t i = 0;
    register const char *s = str;
    register char c;
    for (c = *s; c != '\0'; c = *++s) {
        i *= 4;
        i += sm::code[c] - '0';
    }
    return i;
}

// Low-quality phred score counter. Returns number of bases with a quality
// score below 20.
int lq_count(const char *str, int len)
{
    int lq = 0;
    for (int i = 0; i < len; i++) {
        int phred = str[i] - 33;
        if (phred < 20)
            lq++;
    }
    return lq;
}

// In-place conversion of sequence to its reverse complement.
void revcomp(char seq[], int len)
{
    int c, i, j;
    for (i = 0, j = len - 1; i < j; i++, j--) {
        c = seq[i];
        seq[i] = sm::comp[seq[j]];
        seq[j] = sm::comp[c];
    }
}

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
