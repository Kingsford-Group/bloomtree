#ifndef BLOOMTREE_H
#define BLOOMTREE_H

#include <string>
#include <queue>
#include "Heap.h"
#include "BF.h"

// this is the max number of BF allowed in memory at once.
extern int BF_INMEM_LIMIT;

class BloomTree {
public:
    BloomTree(const std::string & f, HashPair hp, int nh);
    ~BloomTree();
    std::string name() const;

    BloomTree* child(int which) const;
    void set_child(int which, BloomTree* c);

    BF* bf() const;

    //BloomTree* union_bloom_filters(const std::string & new_name, BloomTree* f2) const;

    int usage() const;
    void increment_usage();

private:
    bool load() const;
    void unload() const;

    static Heap<const BloomTree> bf_cache;

    std::string filename;
    HashPair hashes;
    int num_hash;
    mutable BF* bloom_filter;
    mutable Heap<const BloomTree>::heap_reference heap_ref;

    BloomTree* children[2];
    BloomTree* parent;
    int usage_count;
};

HashPair* get_hash_function(const std::string & matrix_file, int & nh);
BloomTree* read_bloom_tree(const std::string & filename);
void write_bloom_tree(const std::string & outfile, BloomTree* root, const std::string & matrix_file);

#endif
