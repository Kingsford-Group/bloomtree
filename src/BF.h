#ifndef BF_H
#define BF_H

#include <string>
#include <sys/mman.h>
#include <sdsl/bit_vectors.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include "Kmers.h"

using HashPair = jellyfish::hash_pair<jellyfish::mer_dna>;

// a kmer bloom filter
class BF {
public:
    BF(const std::string & filename, HashPair hp, int nh);
    ~BF();

    void load();

    bool contains(const jellyfish::mer_dna & m) const;
    bool contains(const std::string & str) const;

    BF* union_with(const std::string & new_name, const BF* f2) const;

private:
    std::string filename;
    sdsl::rrr_vector<255>* bits;

    HashPair hashes;
    unsigned long num_hash;
};


#endif
