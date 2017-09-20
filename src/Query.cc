#include "Query.h"
#include "Kmers.h"
#include "util.h"
#include <cassert>

float QUERY_THRESHOLD = 0.9;

// ** THIS IS NOW PARTIALLY DEPRICATED. ONLY WORKS WITH HARDCODED SIMILARITY TYPE
void assert_is_union(BloomTree* u) {
    BloomTree::protected_cache(true);
    BF* ubf = u->child(0)->bf()->union_with("union", u->child(1)->bf());
    BloomTree::protected_cache(false);

    // if the similaritity isn't 100%, we stop
    std::ostringstream oss;
    uint64_t sim = ubf->similarity(u->bf(),0);
    if (sim != ubf->size()) {
        std::cerr << "Filter at " << u->name() << " is not the union of its two children!" << std::endl;
        std::cerr << "Sim= " << sim << "Size= " << ubf->size() << std::endl;
        std::cerr << "Children:" << u->child(0)->name() << " " << u->child(1)->name() << std::endl;
        std::cerr << std::endl;
    } else {
        std::cerr << "Filter at " << u->name() << " looks good." << std::endl;
        std::cerr << "Children:" << u->child(0)->name() << " " << u->child(1)->name() << std::endl;
        std::cerr << std::endl;
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


void draw_bt_recur(BloomTree* root, std::ostream& out) {
    if (root == nullptr) return;

    std::string current = quote(test_basename(root->name(), ".bf.bv"));
    
    if (root->child(0)) {
        out << current << " -> " << quote(test_basename(root->child(0)->name(), ".bf.bv")) << " ; " << std::endl;
        draw_bt_recur(root->child(0), out);
    } 
    if (root->child(1)) {
        out << current << " -> " << quote(test_basename(root->child(1)->name(), ".bf.bv")) << " ; " << std::endl;
        draw_bt_recur(root->child(1), out);
    } 
}

void draw_bt(BloomTree* root, std::string outfile) {
    std::ofstream out(outfile.c_str());
    DIE_IF(!out, "Couldn't open output file");
    out << "digraph BloomTree {" << std::endl;
    draw_bt_recur(root, out);
    out << "}" << std::endl;
}

void compress_bt(BloomTree* root) {
	if (root == nullptr) return;

	root->bf()->compress();

	if (root->child(0)) {
		compress_bt(root->child(0));
	}
	if (root->child(1)){
		compress_bt(root->child(1));
	}
}

bool query_passes(BloomTree* root, const std::set<jellyfish::mer_dna> & q) {
    auto bf = root->bf();
    unsigned c = 0;
    for (const auto & m : q) {
        //DEBUG: std::cout << "checking: " << m.to_str();
        if (bf->contains(m)) c++;
        //DEBUG: std::cout << c << std::endl;
    }
    return (c >= QUERY_THRESHOLD * q.size());
}


// return true if the filter at this node contains > QUERY_THRESHOLD kmers
bool query_passes(BloomTree* root, QueryInfo*  q) {//const std::set<jellyfish::mer_dna> & q) {
    float weight = 1.0;
    auto bf = root->bf();
    float c = 0;
    unsigned n = 0;
    bool weighted = 0;
    if (q->weight.empty()){
	weighted=0;
    } else { weighted = 1; }
    for (const auto & m : q->query_kmers) {
        //DEBUG: std::cout << "checking: " << m.to_str();
	if (weighted){
		if(q->weight.size() > n){ 
			weight=q->weight[n]; 
		} else {
			std::cerr << "Number of weights (" << q->weight.size() <<") less than query kmers (" << q->query_kmers.size() << ")."  << std::endl;
			exit(3);
		}
	}
        if (bf->contains(m)) c+=weight;
	n++;
        //DEBUG: std::cout << c << std::endl;
    }
    std::cerr << root->name() <<std::endl;
    std::cerr << c << " " << QUERY_THRESHOLD << " " << q->query_kmers.size() << std::endl;
    return (c >= QUERY_THRESHOLD * q->query_kmers.size());
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
        if (query_passes(root, q)) { //q->query_kmers)) {
            if (has_children) {
                pass.emplace_back(q);
            } else {
                q->matching.emplace_back(root);
                n++;
            }
        } 
    }

    // $(node name) $(internal / leaf) $(number of matches)
    if (has_children) { //Changing format
        std::cout << root->name() << " internal " << pass.size() << std::endl;
    } else {
        std::cout << root->name() << " leaf " << n << std::endl;
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


void query_leaves (BloomTree* root, QuerySet & qs) {
    // how many children do we have?
    bool has_children = root->child(0) || root->child(1);

    // construct the set of queries that pass this node
	// But only for leaf nodes
	unsigned n=0;
	if (!has_children) {
    		for (auto & q : qs) {
		        if (query_passes(root, q)) {
		                q->matching.emplace_back(root);
		                n++;
       	 		}
    		}
	}

    // $(node name) $(internal / leaf) $(number of matches)
    if (!has_children) { //Changing format	
        std::cout << root->name() << " leaf " << n << std::endl;
    }

        // if present, recurse into left child
        if (root->child(0)) {
            query_leaves(root->child(0), qs);
        }

        // if present, recurse into right child
        if (root->child(1)) {
            query_leaves(root->child(1), qs);
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

void batch_weightedquery_from_file(
    BloomTree* root,
    const std::string & fn,
    const std::string & wf,
    std::ostream & o
) {
    // read in the query lines from the file.
    std::string line;
    std::string wfline; 
    QuerySet qs;
    std::ifstream in(fn);
    std::ifstream wfin(wf);
    DIE_IF(!in.good(), "Couldn't open query file.");
    DIE_IF(!wfin.good(), "Couldn't open weight file.");
    std::size_t n = 0;
    while (getline(in, line)) {
	getline(wfin, wfline);
        line = Trim(line);
	wfline = Trim(wfline);
	
        if (line.size() < jellyfish::mer_dna::k()) continue;
        qs.emplace_back(new QueryInfo(line, wfline));
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


void leaf_query_from_file(
	BloomTree* root,
	const std::string & fn,
	std::ostream & o
) {
	std::string line;
	QuerySet qs;
	std::ifstream in(fn);
	DIE_IF(!in.good(), "Couldn't open query file.");
	std::size_t n=0;
	while (getline(in, line)) {
		line = Trim(line);
		if (line.size() < jellyfish::mer_dna::k()) continue;
		qs.emplace_back(new QueryInfo(line));
		n++;
	}
	in.close();
	std::cerr << "Read " << n << " queries." << std::endl;

	// batch process the queries on ONLY the leaves
	query_leaves(root, qs);
	print_query_results(qs, o);
	
	for (auto & p : qs) {
		delete p;
	}
}
