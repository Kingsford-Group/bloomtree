#ifndef COUNT_H
#define COUNT_H

#include "BF.h"
#include <string>

bool count(
    std::string infilen,
    std::string outfilen,
    HashPair hp,
    int nh,
    int num_threads = 16
    );
#endif
