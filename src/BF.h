#ifndef BF_H
#define BF_H

#include <string>
#include <sys/mman.h>
#include <sdsl/bit_vectors.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include "Kmers.h"

typedef jellyfish::hash_pair<jellyfish::mer_dna> HashPair;

// a kmer bloom filter
class BF {
public:
    BF(const std::string & filename, HashPair hp);
    ~BF();

    void load();
    void save();

    bool contains(const jellyfish::mer_dna & m);
    bool contains(const std::string & str);

private:
    std::string filename;
    std::string matrix_file;
    sdsl::rrr_vector<127>* bits;

    unsigned long num_hash;
    HashPair hashes;
};


#endif
