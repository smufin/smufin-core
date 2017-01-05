#ifndef __SM_UTIL_H__
#define __SM_UTIL_H__

#include <functional>
#include <vector>

template<typename T>
void spawn(std::string name, std::function<void(T)> func, std::vector<T> list);
void spawn(std::string name, std::function<void(int)> func, int n);

float estimate_sparse(uint64_t n, size_t k, size_t v);

#endif
