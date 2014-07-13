#ifndef QUERY_H
#define QUERY_H

#include <set>
#include <string>
#include <vector>
#include <iostream>

#include "BloomTree.h"
void query_from_file(
    BloomTree* root, 
    const std::string & fn,
    std::ostream & o
);

void query_string(
    BloomTree* root, 
    const std::string & q,
    std::vector<BloomTree*> & out
);

void query(
    BloomTree* root, 
    const std::set<jellyfish::mer_dna> & q, 
    std::vector<BloomTree*> & out
);

#endif
