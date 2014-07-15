#include "BF.h"
#include "Kmers.h"
#include "util.h"

#include <jellyfish/file_header.hpp>

BF::BF(const std::string & f, HashPair hp, int nh) :
    filename(f),
    bits(nullptr),
    hashes(hp),
    num_hash(nh)
{ 
    load();
}

BF::~BF() {
    delete bits;
}

// returns true iff the bloom filter contains the given kmer
bool BF::contains(const jellyfish::mer_dna & m) const {
    uint64_t h0 = hashes.m1.times(m);
    uint64_t h1 = hashes.m2.times(m);

    //DEBUG: std::cout << "size = " << bits->size() << std::endl;
    
    const size_t base = h0 % bits->size();
    const size_t inc = h1 % bits->size();

    for (unsigned long i = 0; i < num_hash; ++i) {
        const size_t pos = (base + i * inc) % bits->size();
        //DEBUG: std::cout << "pos=" << pos << std::endl;
        if ((*bits)[pos] == 0) return false;
    }
    return true;
}

// convience function
bool BF::contains(const std::string & str) const {
    return contains(jellyfish::mer_dna(str));
}


// read the bit vector and the matrices for the hash functions.
void BF::load() {
    // read the actual bits
    bits = new sdsl::rrr_vector<255>();
    sdsl::load_from_file(*bits, filename);
}


// create a new RRR bloom filter that is the union of this BF and the given BF.
// Will re-use the hashes from this and both BFs must use exactly the same hash
// (not checked).
BF* BF::union_with(const std::string & new_name, const BF* f2) const {
    assert(bits->size() == f2->size());

    // create an uncompressed version of this bf
    sdsl::bit_vector b(bits->size(), 0);

    std::cerr << "Performing OR..." << std::endl;

    // union it with f2
    for (unsigned long i = 0; i < b.size(); i++) {
        b[i] = (*bits)[i] | (*f2->bits)[i];
    }

    // create a new BF wraper for the new BF
    BF* out = new BF(new_name, hashes, num_hash);

    // "load" the new BF by converting it to a RRR
    out->bits = new sdsl::rrr_vector<255>(b);
    return out;
}

