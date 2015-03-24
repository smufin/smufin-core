#ifndef __SM_HASH_H__
#define __SM_HASH_H__

#include <stdint.h>
#include <cstring>

uint64_t murmur_hash(const void * key, int len, unsigned int seed);

template<typename T> struct sm_hasher {
    size_t operator()(const T& t) const {
        return murmur_hash(&t, sizeof(t), 0);
    }
};

#endif
