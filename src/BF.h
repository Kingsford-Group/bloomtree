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
    virtual ~BF();

    virtual void load();
    virtual void save();

    virtual int operator[](uint64_t pos) const;
    virtual void set_bit(uint64_t p);
    virtual uint64_t size() const;

    virtual bool contains(const jellyfish::mer_dna & m) const;
    bool contains(const std::string & str) const;

    void add(const jellyfish::mer_dna & m);

    virtual uint64_t similarity(const BF* other, int type) const;
    virtual std::tuple<uint64_t, uint64_t> b_similarity(const BF* other) const;
    virtual BF* union_with(const std::string & new_name, const BF* f2) const;
    virtual void union_into(const BF* f2);
    virtual uint64_t count_ones() const;
    virtual void compress();
protected:
    std::string filename;
    sdsl::rrr_vector<255>* bits;

    HashPair hashes;
    unsigned long num_hash;
};

class UncompressedBF : public BF {
public:
    UncompressedBF(const std::string & filename, HashPair hp, int nh);
    virtual ~UncompressedBF();

    virtual void load();
    virtual void save();

    virtual int operator[](uint64_t pos) const;
    virtual void set_bit(uint64_t p);
    virtual uint64_t size() const;
    virtual uint64_t similarity(const BF* other, int type) const;
    virtual std::tuple<uint64_t, uint64_t> b_similarity(const BF* other) const;
    virtual BF* union_with(const std::string & new_name, const BF* f2) const;
    virtual void union_into(const BF* f2);
    virtual uint64_t count_ones() const;
    virtual void compress();
protected:
    sdsl::bit_vector* bv;
};

sdsl::bit_vector* union_bv_fast(const sdsl::bit_vector & b1, const sdsl::bit_vector& b2);
BF* load_bf_from_file(const std::string & fn, HashPair hp, int nh);

#endif
