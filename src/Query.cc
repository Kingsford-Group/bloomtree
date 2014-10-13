#include "Query.h"
#include "Kmers.h"
#include "util.h"

#include <cassert>

float QUERY_THRESHOLD = 0.9;

// return true if the filter at this node contains > QUERY_THRESHOLD kmers
bool query_passes(BloomTree* root, const std::set<jellyfish::mer_dna> & q) {
    assert(q.size() > 0);
    auto bf = root->bf();
    unsigned c = 0;
    for (const auto & m : q) {
        //DEBUG: std::cout << "checking: " << m.to_str(); 
        if (bf->contains(m)) c++;
        //DEBUG: std::cout << c << std::endl;
    }
    return (c >= QUERY_THRESHOLD * q.size());
}

void assert_is_union(BloomTree* u) {
    BF* ubf = u->child(0)->bf()->union_with("union", u->child(1)->bf());

    // if the similaritity isn't 100%, we stop
    std::ostringstream oss;
    uint64_t sim = ubf->similarity(u->bf());
    if (sim != ubf->size()) {
        std::cerr << "Filter at " << u->name() << " is not the union of its two children!" << std::endl;
        std::cerr << "Sim= " << sim << "Size= " << ubf->size() << std::endl;
        DIE("Stopping.");
    } else {
        std::cerr << "Filter at " << u->name() << " looks good." << std::endl;
    }
    delete ubf;
}

void check_bt(BloomTree* root) {
    if (root == nullptr) return;

    if (root->child(0) && root->child(1)) {
        assert_is_union(root);
    }

    check_bt(root->child(0));
    check_bt(root->child(1));
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
    DIE_IF(!in.good(), "Couldn't open query file.");
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

/******
 * Batch querying
 ******/

void print_query_results(const QuerySet & qs, std::ostream & out) {
    for (auto& q : qs) {
        out << "*" << q->query << " " << q->matching.size() << std::endl;
        for (const auto& n : q->matching) {
            out << n->name() << std::endl;
        }
    }
}


void query_batch(BloomTree* root, QuerySet & qs) {
    // how many children do we have?
    bool has_children = root->child(0) || root->child(1);

    // construct the set of queries that pass this node
    QuerySet pass;
    unsigned n = 0;
    for (auto & q : qs) {
        if (query_passes(root, q->query_kmers)) {
            if (has_children) {
                pass.emplace_back(q);
            } else {
                q->matching.emplace_back(root);
                n++;
            }
        } 
    }

    if (has_children) {
        std::cerr << "At " << root->name() << ", " << pass.size() << " queries passed." << std::endl;
    } else {
        std::cerr << "At leaf " << root->name() << ", " << n << " queries matched." << std::endl;
    }

    if (pass.size() > 0) {
        // if present, recurse into left child
        if (root->child(0)) {
            query_batch(root->child(0), pass);
        }

        // if present, recurse into right child
        if (root->child(1)) {
            query_batch(root->child(1), pass);
        }
    }
} 


void batch_query_from_file(
    BloomTree* root, 
    const std::string & fn,
    std::ostream & o
) { 
    // read in the query lines from the file.
    std::string line;
    QuerySet qs;
    std::ifstream in(fn);
    DIE_IF(!in.good(), "Couldn't open query file.");
    std::size_t n = 0;
    while (getline(in, line)) {
        line = Trim(line);
        if (line.size() < jellyfish::mer_dna::k()) continue;
        qs.emplace_back(new QueryInfo(line));
        n++;
    }
    in.close();
    std::cerr << "Read " << n << " queries." << std::endl;

    // batch process the queries
    query_batch(root, qs);
    print_query_results(qs, o);

    // free the query info objects
    for (auto & p : qs) {
        delete p;
    }
}

