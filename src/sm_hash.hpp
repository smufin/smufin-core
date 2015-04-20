#ifndef __SM_HASH_H__
#define __SM_HASH_H__

#include <stdint.h>
#include <string>
#include <cstring>

uint64_t murmur_hash(const void * key, int len, unsigned int seed);

template<typename T> struct sm_hasher {
    size_t operator()(const T& t) const {
        return murmur_hash(&t, sizeof(t), 0);
    }
};

template<> struct sm_hasher<std::string> {
    size_t operator()(const std::string& t) const {
        return murmur_hash(t.c_str(), t.size(), 0);
    }
};

struct eqstr {
    bool operator()(const char* s1, const char* s2) const {
        return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
    }
};

#endif
