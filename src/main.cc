#include "Query.h"
#include "BloomTree.h"
#include "BF.h"
#include "util.h"

#include <string>
#include <getopt.h>
#include <jellyfish/file_header.hpp>

/* TODO
 */

// various commandline filenames
std::string command;
std::string matrix_file;
std::string bloom_tree_file;
std::string query_file;
std::string out_file;
std::string jfbloom_file;

const char * OPTIONS = "";

static struct option LONG_OPTIONS[] = {
     {0,0,0,0}
};

void print_usage() {
    std::cerr 
        << "Usage: bt [query|convert] ...\n"
        << "    \"query\" matrixfile bloomtreefile queryfile\n"
        << "    \"convert\" jfbloomfilter outfile\n"
        << std::endl;
    exit(3);
}

int process_options(int argc, char* argv[]) {
    int a;
    while ((a=getopt_long(argc, argv, OPTIONS, LONG_OPTIONS, 0)) != -1) {
        switch(a) {
            default:
                std::cerr << "Unknown option." << std::endl;
                print_usage();
        }
    }

    if (optind >= argc) print_usage();
    command = argv[optind];
    if (command == "query") {
        if (optind >= argc-3) print_usage();
        matrix_file = argv[optind+1];
        bloom_tree_file = argv[optind+2];
        query_file = argv[optind+3];

    } else if (command == "convert") {
        if (optind >= argc-2) print_usage();
        jfbloom_file = argv[optind+1];
        out_file = argv[optind+2];
    }
    return optind;
}

// non-zero if bit is set
inline char bit(char * buf, unsigned long bit) {
    char byte = buf[bit / 8];
    char bit_mask = 1 << (bit % 8); 
    return byte & bit_mask;
}

int main(int argc, char* argv[]) {
    std::cerr << "Starting Bloom Tree" << std::endl;

    process_options(argc, argv);

    if (command == "query") {
        std::cerr << "Loading hashes from " << matrix_file << std::endl; 
        std::ifstream in(matrix_file.c_str(), std::ios::in | std::ios::binary);
        jellyfish::file_header header(in);
        DIE_IF(!in.good(), "Couldn't parse bloom filter header!");
        HashPair hashes(header.matrix(1), header.matrix(2));
        in.close();

        auto k = header.key_len();
        jellyfish::mer_dna::k(k / 2);
        std::cerr << "Read hashes for k=" << jellyfish::mer_dna::k() 
            << std::endl;
        std::cerr << "# Hash applications=" << header.nb_hashes() << std::endl;

        std::cerr << "Loading bloom tree topology: " << bloom_tree_file 
            << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file, hashes, header.nb_hashes());

        std::cerr << "Querying..." << std::endl;
        query_from_file(root, query_file, std::cout);

    } else if (command == "convert") {
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
        
        // covert that raw bit vector into an rrr compressed vector &
        // save it.
        sdsl::rrr_vector<255> rrr(b);
        std::cerr << "Compressed RRR vector is " 
            << sdsl::size_in_mega_bytes(rrr) << std::endl;
        sdsl::store_to_file(rrr, out_file);

        delete buf;
    }
    std::cerr << "Done." << std::endl;
}
