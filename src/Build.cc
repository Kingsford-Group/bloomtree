#include "Build.h"
#include "BloomTree.h"
#include "util.h"

#include <sys/mman.h>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/file_header.hpp>
#include <sdsl/bit_vectors.hpp>

// read a list of bloom filters, 1 per line from the given file
std::vector<std::string> read_filter_list(const std::string & inf) {
    std::ifstream in(inf);
    std::vector<std::string> v;
    std::string line;
    while(getline(in, line)) {
        v.emplace_back(Trim(line));
    }
    return v;
}

// compute the index of the children of position i
std::size_t complete_tree_child(std::size_t i, unsigned child) {
    switch (child) {
        case 0: return 2*i + 1;
        case 1: return 2*i + 2;
        default: DIE("invalid child in complete_tree_child");
    }
    return 0;
}


// builds a near-complete bt (the last level might be only partially full)
// with the given bf as the leaves
void build_bloom_tree_filters(
    const std::vector<std::string> & fs, 
    const std::string & matrix_file,
    const std::string & outf
) {
    // create the hashes
    int nh = 0;
    HashPair* hashes = get_hash_function(matrix_file, nh); 

    // compute the number of nodes needed:
    unsigned lp2 = int(log2(fs.size()));  // largest power of 2 <= n
    unsigned nb_nodes = (2*lp2 - 1) + 2*(fs.size()-lp2);

    std::cerr << "Tree will have " << nb_nodes << " nodes" << std::endl;
    // v holds the nodes of the semi-complete tree
    std::vector<BloomTree*> v(nb_nodes, nullptr);

    // create the leaves
    std::cerr << "Creating leaves..." << std::endl;
    int i = int(v.size()-1); // must be signed for subsequent loop!
    for (const auto & s : fs) {
        v[i--] = new BloomTree(s, *hashes, nh);
    }

    // i points to the first node we must union
    char buf[100]; // c++ fail here
    for (; i >= 0; i--) {
        sprintf(buf, "union%d.rrr", i);
        std::string new_name = buf;
        std::cerr << "Creating union: " << i << " " << new_name << std::endl;

        v[i] = v[complete_tree_child(i, 0)]->union_bloom_filters(
            new_name,
            v[complete_tree_child(i, 1)]
        );
    }

    write_bloom_tree(outf, v[0], matrix_file);

    // free the memory for the bloom nodes
    for (auto & p : v) {
        delete p;
    }
}

// non-zero if bit is set
inline char bit(char * buf, unsigned long bit) {
    char byte = buf[bit / 8];
    char bit_mask = 1 << (bit % 8); 
    return byte & bit_mask;
}

void convert_jfbloom_to_rrr(const std::string & jfbloom_file, const std::string & out_file) {
    // Load the JF Bloom Filter
    std::cerr << "Reading: " << jfbloom_file << std::endl;
    std::ifstream in(jfbloom_file.c_str(), std::ios::in|std::ios::binary);
    DIE_IF(!in.good(), "Couldn't open jellyfish bloom filter");
    jellyfish::file_header header(in);

    // length = # of bits; num_bytes = # of bytes
    auto length = header.size();
    auto num_bytes = length / 8 + (length % 8 != 0);

    std::cerr << "BF length = " << length << " bits, occupying "
        << num_bytes << " bytes" << std::endl; 

    // suck in the bits from the bf
    char * buf = new char[num_bytes];
    in.read(buf, num_bytes);
    in.close();
    
    // Create a new SDSL bitvector
    sdsl::bit_vector b(length, 0);
    for (unsigned long i = 0; i < length; i++) {
        if (bit(buf, i)) {
            b[i] = 1;
        }
    }
    
    // covert that raw bit vector into an rrr compressed vector & save it.
    sdsl::rrr_vector<255> rrr(b); std::cerr << "Compressed RRR vector is " <<
        sdsl::size_in_mega_bytes(rrr) << std::endl; sdsl::store_to_file(rrr,
                out_file);

    delete buf;
}
