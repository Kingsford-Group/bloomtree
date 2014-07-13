#ifndef KMERS_H
#define KMERS_H
#include <set>
#include <string>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/rectangular_binary_matrix.hpp>

//int acgt(char c);
//Kmer kmer_to_bits(const std::string & str);
std::set<jellyfish::mer_dna> kmers_in_string(const std::string & str);

#endif
