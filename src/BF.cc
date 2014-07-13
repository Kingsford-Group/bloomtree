#include "BF.h"
#include "Kmers.h"
#include "util.h"

#include <jellyfish/file_header.hpp>

BF::BF(const std::string & f, HashPair hp, int nh) :
    filename(f),
    bits(0),
    hashes(hp),
    num_hash(nh)
{ 
    load();
}

BF::~BF() {
    delete bits;
}

bool BF::contains(const jellyfish::mer_dna & m) {
    uint64_t h0 = hashes.m1.times(m);
    uint64_t h1 = hashes.m2.times(m);
    
    const size_t base = h0 % bits->size();
    const size_t inc = h1 % bits->size();

    for (unsigned long i = 0; i < num_hash; ++i) {
        const size_t pos = (base + i * inc) % bits->size();
        if ((*bits)[pos] == 0) return false;
    }
    return true;
}

// convience function
bool BF::contains(const std::string & str) {
    return contains(jellyfish::mer_dna(str));
}

// save the bitvector to a file
void BF::save() {
    sdsl::store_to_file(*bits, filename); 
}

// read the bit vector and the matrices for the hash functions.
void BF::load() {
    // read the actual bits
    sdsl::load_from_file(*bits, filename);
}



