#include "Count.h"
#include "Kmers.h"
#include "BF.h"
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/file_header.hpp>
#include <sdsl/bit_vectors.hpp>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>

#include <vector>

/*==== COPIED FROM THE JF count_dump example ====*/

typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna> mer_hash_type;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<std::vector<std::string>::iterator > > sequence_parser_type;
typedef jellyfish::mer_iterator<sequence_parser_type, jellyfish::mer_dna> mer_iterator_type;


class mer_counter : public jellyfish::thread_exec {
  mer_hash_type&                    mer_hash_;
  jellyfish::stream_manager<std::vector<std::string>::iterator > streams_;
  sequence_parser_type              parser_;
  const bool                        canonical_;

public:
  mer_counter(int nb_threads, mer_hash_type& mer_hash,
              std::vector<std::string>::iterator file_begin, std::vector<std::string>::iterator file_end,
              bool canonical)
    : mer_hash_(mer_hash)
    , streams_(file_begin, file_end)
    , parser_(jellyfish::mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
    , canonical_(canonical)
  { }

  virtual void start(int thid) {
    mer_iterator_type mers(parser_, canonical_);

    for( ; mers; ++mers)
      mer_hash_.add(*mers, 1);
    mer_hash_.done();
  }
};

/*=== END COPY ===*/


// run JF count, and build a BF and save it to disk

enum OPERATION { COUNT, PRIME, UPDATE };

bool count(
    std::string infilen,
    std::string outfilen,
    HashPair hp,
    int nh,
    uint64_t bf_size,
    int num_threads,
    unsigned cutoff_count
    ) {
    
    // jellyfish default counting values
    const uint64_t hash_size = 10000000;
    const uint32_t num_reprobes = 126;
    const uint32_t counter_len = 7;
    const bool canonical = true;

    // create the hash
    mer_hash_type mer_hash(hash_size, jellyfish::mer_dna::k()*2, counter_len, num_threads, num_reprobes);

    // create a mock up of the array of file names
    std::vector<std::string> files;
    files.push_back(infilen);

    // count the kmers
    mer_counter counter(num_threads, mer_hash, files.begin(), files.end(), canonical);
    counter.exec_join(num_threads);

    // build the BF
    UncompressedBF bf(outfilen, hp, nh, bf_size);

    // add each kmer to the BF
    const auto jf_ary = mer_hash.ary();
    const auto end = jf_ary->end();
    std::cerr << "Right before cutoff count: " << cutoff_count << std::endl;
    for(auto kmer = jf_ary->begin(); kmer != end; ++kmer) {
        auto& key_val = *kmer;
        if (key_val.second >= cutoff_count) {
            bf.add(key_val.first);
        }
    }
    bf.save();
    return true;
}


// bt count k FASTA OUTFILE

 /*
    // create the hash
    jellyfish::mer_hash ary(hash_size, k*2, counter_len, num_threads, num_reprobes);

    // prime the kmer values
    jellyfish::mer_counter primer(num_threads, ary,
            files.begin(), files.end(),
            files.end(), files.end(), 1, PRIME); 
    primer.exec_join(num_threads);

    // count the kmers
    jellyfish::mer_counter counter(num_threads, ary,
            files.begin(), files.end(),
            files.end(), files.end(),
            1, UPDATE);
    counter.exec_join(num_threads);
    */
