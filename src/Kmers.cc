#include "Kmers.h"
#include "util.h"

/*
// return the number for each DNA base
int acgt(char c) {
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default:
        DIE("Unknown nucleotide");
        return 0;
    }
}

// convert a kmer to a uint64_t
Kmer kmer_to_bits(const std::string & str) {
    Kmer b;
    for (int i = 0; i < str.size(); i++) {
        b = (b<<2) | acgt(str[i]);
    }
    return b;
}
*/

std::set<jellyfish::mer_dna> kmers_in_string(const std::string & str) {
    auto k = jellyfish::mer_dna::k();
    set<jellyfish::mer_dna> s;
    for (int i = 0; i <= str.size() - k; i++) {
        s.insert(jellyfish::mer_dna(str.substr(i, k)));
    }
    return s;
}
