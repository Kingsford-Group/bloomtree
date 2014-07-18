#include "Query.h"
#include "Kmers.h"
#include "util.h"

float QUERY_THRESHOLD = 0.8;

// return true if the filter at this node contains > QUERY_THRESHOLD kmers
bool query_passes(BloomTree* root, const std::set<jellyfish::mer_dna> & q) {
    auto bf = root->bf();
    int c = 0;
    for (const auto & m : q) {
        //DEBUG: std::cout << "checking: " << m.to_str(); 
        if (bf->contains(m)) c++;
        //DEBUG: std::cout << c << std::endl;
    }
    return (c > QUERY_THRESHOLD * q.size());
}

// recursively walk down the tree, proceeding to children only
// if their parent passes the query threshold; 
void query(
    BloomTree* root, 
    const std::set<jellyfish::mer_dna> & q, 
    std::vector<BloomTree*> & out
) {
    root->increment_usage();
    if (query_passes(root, q)) {
        //DEBUG: std::cout << "passed at " << root->name() << std::endl;
        int children = 0;
        if (root->child(0)) {
            query(root->child(0), q, out);
            children++;
        }
        if (root->child(1)) {
            query(root->child(1), q, out);
            children++;
        }
        if (children == 0) {
            out.push_back(root);
        }
    }
}

// same as query() but the string is first converted into a set of kmers.
void query_string(
    BloomTree* root, 
    const std::string & q,
    std::vector<BloomTree*> & out
) {
    query(root, kmers_in_string(q), out);
}

// read 1 query per line, execute it, and print to the output stream o the
// results in the format:
//      *QUERY number_results
//      BF names
//
void query_from_file(
    BloomTree* root, 
    const std::string & fn,
    std::ostream & o
) { 
    std::vector<BloomTree*> out;
    std::string line;

    std::ifstream in(fn);
    while (getline(in, line)) {
        line = Trim(line);
        if (line.size() < jellyfish::mer_dna::k()) continue;

        o << "*" << line;
        query_string(root, line, out);
        o << " " << out.size() << std::endl;

        for (const auto& n : out) {
            o << n->name() << std::endl;
        }

        out.clear();
    }
}
