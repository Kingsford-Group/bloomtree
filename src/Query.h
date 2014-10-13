#ifndef QUERY_H
#define QUERY_H

#include <set>
#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "BloomTree.h"

extern float QUERY_THRESHOLD;

struct QueryInfo {
    QueryInfo(const std::string & q) : query(q), query_kmers(kmers_in_string(q)) {}
    ~QueryInfo() {}

    std::string query;
    std::set<jellyfish::mer_dna> query_kmers;
    std::vector<const BloomTree*> matching;
};

using QuerySet = std::list<QueryInfo*>;

void query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o); 
void query_string(BloomTree* root, const std::string & q, std::vector<BloomTree*> & out);
void query(BloomTree* root, const std::set<jellyfish::mer_dna> & q, std::vector<BloomTree*> & out);
void check_bt(BloomTree* root);

#endif
