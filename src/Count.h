#ifndef COUNT_H
#define COUNT_H

#include "BF.h"
#include <string>

bool count(
    std::string infilen,
    std::string outfilen,
    HashPair hp,
    int nh,
    uint64_t bf_size,
    int num_threads = 16,
    unsigned cutoff_count = 3
    );
#endif
