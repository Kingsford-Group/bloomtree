#ifndef BLOOMTREE_H
#define BLOOMTREE_H

#include <string>
#include <queue>
#include "Heap.h"
#include "BF.h"

// this is the max number of BF allowed in memory at once.
static const int BF_INMEM_LIMIT = 32;

class BloomTree {
public:
    BloomTree(const std::string & f, HashPair hp, int nh);
    ~BloomTree();
    BloomTree* child(int which);
    void set_child(int which, BloomTree* c);
    BF* bf();

    std::string name();

    int usage();
    void increment_usage();

private:
    bool load();
    void unload();

    static Heap<BloomTree> bf_cache;

    std::string filename;
    HashPair hashes;
    int num_hash;
    BF* bloom_filter;
    BloomTree* children[2];
    BloomTree* parent;
    int usage_count;
    Heap<BloomTree>::heap_reference heap_ref;
};

BloomTree* read_bloom_tree(const std::string & filename);

#endif
