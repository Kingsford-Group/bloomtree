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
    int num_children() const;
    void set_parent(const BloomTree* p);
    const BloomTree* get_parent() const;
    uint64_t similarity(BloomTree* other) const;
    std::tuple<uint64_t, uint64_t> b_similarity(BloomTree* other) const;
    BF* bf() const;

    BloomTree* union_bloom_filters(const std::string & new_name, BloomTree* f2);
    void union_into(const BloomTree* other);

    int usage() const;
    void increment_usage() const;
    static void protected_cache(bool b);

private:
    bool load() const;
    void unload() const;

    static Heap<const BloomTree> bf_cache;
    static void drain_cache();

    std::string filename;
    HashPair hashes;
    int num_hash;
    mutable BF* bloom_filter;
    mutable Heap<const BloomTree>::heap_reference* heap_ref;

    BloomTree* children[2];
    BloomTree* parent;
    mutable int usage_count;
    mutable bool dirty;
};

HashPair* get_hash_function(const std::string & matrix_file, int & nh);
BloomTree* read_bloom_tree(const std::string & filename);
void write_bloom_tree(const std::string & outfile, BloomTree* root, const std::string & matrix_file);

#endif
