#include "BloomTree.h"
#include "util.h"
#include "BF.h"
#include "gzstream.h"

#include <fstream>
#include <list>
#include <cassert>
#include <jellyfish/file_header.hpp>

Heap<const BloomTree> BloomTree::bf_cache;
int BF_INMEM_LIMIT = 100;

// construct a bloom filter with the given filter backing.
BloomTree::BloomTree(
    const std::string & f, 
    HashPair hp,
    int nh
) :
    filename(f),
    hashes(hp),
    num_hash(nh),
    bloom_filter(0),
    heap_ref(nullptr),
    parent(0),
    usage_count(0),
    dirty(false)
{
    children[0] = nullptr;
    children[1] = nullptr;
}

// free the memory for this node.
BloomTree::~BloomTree() {
    unload();
}

std::string BloomTree::name() const {
    return filename;
}

// Return the node for the given child
BloomTree* BloomTree::child(int which) const { 
    assert(which >= 0 || which < 2);
    return children[which]; 
}

// Set the given child
void BloomTree::set_child(int which, BloomTree* c) {
    assert(which >= 0 || which < 2);
    c->parent = this;
    children[which] = c;
}

int BloomTree::num_children() const {
    return ((children[0]==nullptr)?0:1) + ((children[1]==nullptr)?0:1);
}

const BloomTree* BloomTree::get_parent() const {
    return this->parent;
}

void BloomTree::set_parent(const BloomTree* p) {
    parent = const_cast<BloomTree*>(p);
}

// return the bloom filter, loading first if necessary
BF* BloomTree::bf() const {
    load();
    return bloom_filter;
}

// return the number of times this bloom filter has been used.
int BloomTree::usage() const {
    return usage_count;
}

// increment the usage counter, and update the filter's location in
// the heap if needed.
void BloomTree::increment_usage() const {
    usage_count++;
    // if we're in the cache, let the cache know we've been used.
    if (heap_ref != nullptr) {
        bf_cache.increase_key(heap_ref, usage_count);
    }
}

// Frees the memory associated with the bloom filter
void BloomTree::unload() const { 
    // you can't unload something until you remove it from the cache
    // DEBUG std::cerr << "Unloading " << name() << std::endl;
    
    // free the memory
    if (bloom_filter != nullptr) {
        if (dirty) {
            bloom_filter->save();
        }
        delete bloom_filter; 
        bloom_filter = nullptr; 
    }
    dirty = false;
}

void BloomTree::drain_cache() {
    // if the cache is too big
    while (bf_cache.size() >= BF_INMEM_LIMIT && !bf_cache.is_protected()) {
        // toss the bloom filter with the lowest usage
        const BloomTree* loser = bf_cache.pop();
        loser->heap_ref = nullptr;

        //std::cerr << "Unloading BF: " << loser->filename   
        //          << " cache size = " << bf_cache.size() << std::endl;
        loser->unload();
    }
}

void BloomTree::protected_cache(bool b) {
    bf_cache.set_protected(b);
    if (!b) {
        BloomTree::drain_cache();
    }
}

// Loads the bloom filtering into memory
bool BloomTree::load() const {
    if (bloom_filter == nullptr) {
        //std::cerr << "Loading BF: " << filename << std::endl;

        // if the cache isn't protected from deleting elements, remove enough
        // elements so that there is 1 cache spot free (if the cache is
        // protected, we're allowed to go over the cache limit)
        if(!bf_cache.is_protected()) BloomTree::drain_cache();

        // read the BF file and set bloom_filter
        bloom_filter = load_bf_from_file(filename, hashes, num_hash);
        bloom_filter->load();
        heap_ref = bf_cache.insert(this, usage());
        dirty = false;

        // since we had to load, we bump up the usage to make it less likely we
        // load again in the near future.
        increment_usage();
    }
    return true;
}


uint64_t BloomTree::similarity(BloomTree* other, int type) const {
    protected_cache(true);
    uint64_t sim = this->bf()->similarity(other->bf(), type);
    protected_cache(false);
    return sim;
}

std::tuple<uint64_t,uint64_t> BloomTree::b_similarity(BloomTree* other) const{
    protected_cache(true);
std::cerr << "Before \n";
    std::tuple<uint64_t, uint64_t> sim = this->bf()->b_similarity(other->bf());
    protected_cache(false);
std::cerr << "After \n";
    return sim;
}


// Create a new node that is the union of the bloom filters
// in two other nodes;
BloomTree* BloomTree::union_bloom_filters(const std::string & new_name, BloomTree* f2) {
    // move the union op into BloomTree?
    BloomTree* bt = new BloomTree(new_name, hashes, num_hash);

    protected_cache(true);
    bt->bloom_filter = bf()->union_with(new_name, f2->bf()); 

    bt->set_child(0, this);
    bt->set_child(1, f2);
    //bf_cache.insert(bt, bt->usage());
    bt->dirty = true;
    bt->unload();

    protected_cache(false);
    return bt; 
}

void BloomTree::union_into(const BloomTree* other) {
    protected_cache(true);
    bf()->union_into(other->bf());
    dirty = true;
    protected_cache(false);
}

/*
BloomTree* create_union_node(BloomTree* T, BloomTree* N) {
    assert(T != nullptr && N != nullptr);
    assert(T->hashes == N->hashes);
    assert(T->num_hashes == N->num_hashes);

    std::string new_name = "union";

    BloomTree* unionNode = new BloomTree(new_name, T->hashes, T->num_hash);
    unionNode->bloom_filter = T->bloom_filter->union_with(new_name, N->bloom_filter);
    unionNode->set_child(0, T);
    unionNode->set_child(1, N);
    return unionNode;
}
*/

HashPair* get_hash_function(const std::string & matrix_file, int & nh) {
    std::cerr << "Loading hashes from " << matrix_file << std::endl; 
    igzstream in(matrix_file.c_str(), std::ios::in | std::ios::binary);
    jellyfish::file_header header(in);
    DIE_IF(!in.good(), "Couldn't parse bloom filter header!");
    HashPair * hp = new HashPair(header.matrix(1), header.matrix(2));
    in.close();

    nh = header.nb_hashes();
    std::cerr << "# Hash applications=" << nh << std::endl;

    jellyfish::mer_dna::k(header.key_len() / 2);
    std::cerr << "Read hashes for k=" << jellyfish::mer_dna::k() << std::endl;
    return hp;
}


/* read a file that defines the bloom tree structure. The
   file has lines of the form:
    Root
    *Child1
    ***Child3
    ****Child4
    *Child2
   where the '*' indicates the level and where "Child1" etc
   are comma-separated lists of attributes, the first one 
   of which must be the BF filename.

   This function will return a pointer to the root of the
   constructed bloom tree.
*/
BloomTree* read_bloom_tree(const std::string & filename, bool read_hashes) {
    std::ifstream in(filename.c_str());

    std::list<BloomTree*> path;
    BloomTree* tree_root = 0;
    int n = 0;
    // if read_hashes is false, you must promise never to access the bloom filters
    HashPair* hashes = new HashPair; // useless hashpair used if read_hashes is false
    int num_hashes = 0;

    std::string node_info;
    while (getline(in, node_info)) {
        node_info = Trim(node_info);
        if (node_info.size() == 0) continue;
        size_t level = node_info.find_first_not_of("*");
        node_info.erase(0, level);

        // each node info is a comma separated list
        std::vector<std::string> fields;
        SplitString(node_info, ',', fields);
        std::string bf_filename = fields[0];
        //std::cerr << "Reading BN info: " << bf_filename << " level = " << level << std::endl;

        n++;

        BloomTree* bn = nullptr;

        // if we're at the root
        if (path.size() == 0) {
            DIE_IF(level != 0, "Root must start in column 0");
            DIE_IF(tree_root != 0, "Can't set root twice!");

            // set the hash function up
            if (read_hashes) {
                DIE_IF(fields.size() < 2, "Must specify hash file for root.");
                hashes = get_hash_function(fields[1], num_hashes);
            }

            // create the root node
            bn = new BloomTree(bf_filename, *hashes, num_hashes); 
            tree_root = bn;
            
        // if we're adding a child
        } else {
            bn = new BloomTree(bf_filename, *hashes, num_hashes); 

            while (path.size() > level) {
                path.pop_back();
            }
            DIE_IF(level != path.size(), 
                "Must increase level by <= 1");

            if (path.back()->child(0) == nullptr) {
                path.back()->set_child(0, bn);
            } else if (path.back()->child(1) == nullptr) {
                path.back()->set_child(1, bn);
            } else {
                DIE("Tried to add >= 2 children to a node.");
            }
        }
        path.push_back(bn);
    }
    delete hashes;

    std::cerr << "Read " << n << " nodes in Bloom Tree" << std::endl;
    
    return tree_root;
}

void write_bloom_tree_helper(std::ostream & out, BloomTree* root, int level=1) {
    std::string lstr(level, '*');

    for (int i = 0; i < 2; i++) {
        if (root->child(i) != nullptr) {
            out << lstr << root->child(i)->name() << std::endl;
            write_bloom_tree_helper(out, root->child(i), level+1);
        }
    }
}

// write the bloom tree file format in a way that can be read by
// read_bloom_tree()
void write_bloom_tree(
    const std::string & outfile, 
    BloomTree* root, 
    const std::string & matrix_file
) {
    std::cerr << "Writing to " << outfile << std::endl;
    std::ofstream out(outfile.c_str());
    out << root->name() << "," << matrix_file << std::endl;
    write_bloom_tree_helper(out, root);
    std::cerr << "Done." << std::endl;
}

void write_compressed_bloom_tree_helper(std::ostream & out, BloomTree* root, int level=1) {
    std::string lstr(level, '*');

    for (int i = 0; i < 2; i++) {
        if (root->child(i) != nullptr) {
            out << lstr << root->child(i)->name() << ".rrr" << std::endl;
            write_compressed_bloom_tree_helper(out, root->child(i), level+1);
        }
    }
}

// write the bloom tree file format in a way that can be read by
// read_bloom_tree()
void write_compressed_bloom_tree(
    const std::string & outfile,
    BloomTree* root,
    const std::string & matrix_file
) {
    std::cerr << "Writing to " << outfile << std::endl;
    std::ofstream out(outfile.c_str());
    out << root->name() << ".rrr," << matrix_file << std::endl;
    write_compressed_bloom_tree_helper(out, root);
    std::cerr << "Done." << std::endl;
}



