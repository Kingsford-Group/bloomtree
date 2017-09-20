#include "Query.h"
#include "Build.h"
#include "BloomTree.h"
#include "BF.h"
#include "util.h"
#include "Count.h"

#include <string>
#include <cstdlib>
#include <getopt.h>

#include <jellyfish/file_header.hpp>

/* TODO:
 */

// various commandline filenames
std::string command;
std::string bloom_tree_file;
std::string query_file;
std::string out_file;
std::string jfbloom_file;
std::string bvfile1, bvfile2;
int sim_type=0;
std::string bloom_storage;
int leaf_only;
std::string weighted="";
unsigned cutoff_count=3;

std::string hashes_file;
unsigned nb_hashes;
uint64_t bf_size;

unsigned num_threads = 16;
//unsigned parallel_level = 3; // no parallelism by default

const char * OPTIONS = "t:p:f:l:c:w:s:";

static struct option LONG_OPTIONS[] = {
    {"max-filters", required_argument, 0, 'f'},
    {"threads", required_argument, 0, 'p'},
    {"query-threshold", required_argument, 0, 't'},
    {"k", required_argument, 0, 'k'},
    {"sim-type", required_argument, 0, 's'},
    {"leaf-only", required_argument,0,'l'},
    {"cutoff", required_argument,0,'c'},
    {"weighted", required_argument,0,'w'},
    {0,0,0,0}
};

void print_usage() {
    std::cerr 
        << "Usage: bt [query|convert|build] ...\n"
        << "    \"hashes\" [-k 20] hashfile nb_hashes\n"
        << "    \"count\" [--cutoff 3] [--threads 16] hashfile bf_size fasta_in filter_out.bf.bv\n"
        << "    \"build\" [--sim-type 0] hashfile filterlistfile outfile\n"
	    << "    \"compress\" bloomtreefile outfile\n"

        << "    \"check\" bloomtreefile\n"
        << "    \"draw\" bloomtreefile out.dot\n"

        << "    \"query\" [--max-filters 1] [-t 0.8] [-leaf-only 0] [--weighted weightfile] bloomtreefile queryfile outfile\n"

        << "    \"convert\" jfbloomfilter outfile\n"
        << "    \"sim\" [--sim-type 0] bloombase bvfile1 bvfile2\n"
        << std::endl;
    exit(3);
}

// construct a new set of hashes for the current k
// also record k and the number of hash applications
void construct_hashes(std::string & hashesfile, int nh) {
    HashPair hp;
    jellyfish::file_header fh;
    fh.matrix(hp.m1, 1);
    fh.matrix(hp.m2, 2);
    fh.key_len(jellyfish::mer_dna::k() * 2);
    fh.nb_hashes(nh);
    std::ofstream hashesout(hashesfile.c_str());
    fh.write(hashesout);
    hashesout.close();
}


int process_options(int argc, char* argv[]) {
    int k = 20;
    //int sim_type = 0;
    int a;
    while ((a=getopt_long(argc, argv, OPTIONS, LONG_OPTIONS, 0)) != -1) {
        switch(a) {
            case 't':
                QUERY_THRESHOLD = atof(optarg);
                break;
            case 'p':
        		num_threads = unsigned(atoi(optarg));
                //parallel_level = unsigned(atoi(optarg));
                break;
            case 'f':
                BF_INMEM_LIMIT = unsigned(atoi(optarg));
                break;
	        case 'l':
    	    	leaf_only = atoi(optarg);
                break;
            case 'k': 
                k = atoi(optarg);
                break;
    	    case 's':
	        	sim_type = atoi(optarg);
        		break;
    	    case 'c':
    	    	cutoff_count = unsigned(atoi(optarg));
		        break;
	        case 'w':
		        weighted = optarg;
        		break;	
            default:
                std::cerr << "Unknown option." << std::endl;
                print_usage();
        }
    }

    jellyfish::mer_dna::k(k);
    std::cerr << "Kmer size = " << jellyfish::mer_dna::k() << std::endl;


    if (optind >= argc) print_usage();
    command = argv[optind];
    if (command == "query") {
        if (optind >= argc-3) print_usage();
        bloom_tree_file = argv[optind+1];
        query_file = argv[optind+2];
        out_file = argv[optind+3];
        //leaf_only = argv[optind+4];

    } else if (command == "check") {
        if (optind >= argc-1) print_usage();
        bloom_tree_file = argv[optind+1];

    } else if (command == "draw") {
        if (optind >= argc-2) print_usage();
        bloom_tree_file = argv[optind+1];
        out_file = argv[optind+2];

    } else if (command == "sim") {
        if (optind >= argc-3) print_usage();
        jfbloom_file = argv[optind+1];
        bvfile1 = argv[optind+2];
        bvfile2 = argv[optind+3];
	//    sim_type = argv[optind+4];

    } else if (command == "convert") {
        if (optind >= argc-2) print_usage();
        jfbloom_file = argv[optind+1];
        out_file = argv[optind+2];

    } else if (command == "build") {
        if (optind >= argc-3) print_usage();
        hashes_file = argv[optind+1];
        query_file = argv[optind+2];
        out_file = argv[optind+3];
        //bloom_storage = argv[optind+4];
        // sim_type = argv[optind+4];

    } else if (command == "hashes") {
        if (optind >= argc-2) print_usage();
        hashes_file = argv[optind+1];
        nb_hashes = atoi(argv[optind+2]);

    } else if (command == "count") {
        if (optind >= argc-4) print_usage();
        hashes_file = argv[optind+1];
        bf_size = atof(argv[optind+2]);
        query_file = argv[optind+3];
        out_file = argv[optind+4];


    } else if (command == "compress") {
        if (optind >= argc-2) print_usage();
        bloom_tree_file = argv[optind+1];
        out_file = argv[optind+2];
    }
    return optind;
}



int main(int argc, char* argv[]) {
    std::cerr << "Starting Bloom Tree" << std::endl;

    process_options(argc, argv);

    if (command == "query") {
        std::cerr << "Loading bloom tree topology: " << bloom_tree_file 
            << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file);

        std::cerr << "In memory limit = " << BF_INMEM_LIMIT << std::endl;

        std::cerr << "Querying..." << std::endl;
        std::ofstream out(out_file);
	if (leaf_only == 1){
		leaf_query_from_file(root, query_file, out);
	} else if (weighted!="") {
		std::cerr << "Weighted query \n";
		batch_weightedquery_from_file(root, query_file, weighted, out);	
	} else {
	        batch_query_from_file(root, query_file, out);
	}

    } else if (command == "draw") {
        std::cerr << "Drawing tree in " << bloom_tree_file << " to " << out_file << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file, false);
        draw_bt(root, out_file);

    } else if (command == "check") {
        BloomTree* root = read_bloom_tree(bloom_tree_file);
        std::cerr << "Checking tree" << std::endl;
        check_bt(root);

    } else if (command == "sim") {
        // read hash functions
        int num_hash;
        HashPair* hashes = get_hash_function(jfbloom_file, num_hash);

        // read bloom filters
        std::cerr << "Loading BFs:" << bvfile1 << " " << bvfile2 << std::endl;
        BF* bf1 = load_bf_from_file(bvfile1, *hashes, num_hash);
        BF* bf2 = load_bf_from_file(bvfile2, *hashes, num_hash);
        bf1->load();
        bf2->load();

        std::cerr << "Computing Sim..." << std::endl;
        uint64_t test = bf1->similarity(bf2, sim_type);
	//std::tuple<uint64_t, uint64_t> sim = bf1->b_similarity(bf2);
	std::cerr << "Temphead: " << test << std::endl;
	//std::cout << bf1->size() << " " << std::get<0>(sim) << " " << std::get<1>(sim) << std::endl;

	//uint64_t sim = bf1->similarity(bf2);
        //std::cout << bf1->size() << " " << sim << std::endl;
	//std::cout << "Similarity: " << sim << std::endl;
        //std::cout << "Difference: " << bf1->size() - sim << std::endl;
        //std::cout << "Ones: " << bf1->count_ones() << " " << bf2->count_ones() << std::endl;
        //std::cout << "Size: " << bf1->size() << std::endl;

        delete bf1;
        delete bf2;

    } else if (command == "convert") {
        std::cerr << "Converting..." << std::endl;
        convert_jfbloom_to_rrr(jfbloom_file, out_file);

    } else if (command == "hashes") {
        // construct a new hashpair
        construct_hashes(hashes_file, nb_hashes);
    } else if (command == "count") {
        int nh;
        HashPair* hp = get_hash_function(hashes_file, nh);
	std::cerr << "Cutoff Count: " << cutoff_count << std::endl;
        count(query_file, out_file, *hp, nh, bf_size, num_threads, cutoff_count);

    } else if (command == "build") {
        std::cerr << "Building..." << std::endl;
        std::vector<std::string> leaves = read_filter_list(query_file); //not a query file
        //build_bt_from_jfbloom(leaves, out_file, parallel_level);
        dynamic_build(hashes_file, leaves, out_file, sim_type); //std::stoi(sim_type));

    } else if (command == "compress") {
        std::cerr << "Compressing.." << std::endl;
        BloomTree* root = read_bloom_tree(bloom_tree_file, false);
        std::ifstream in(bloom_tree_file.c_str());
        std::string header;
        getline(in, header);
        std::vector<std::string> fields;
        SplitString(header, ',', fields);
            
        compress_bt(root);
        write_compressed_bloom_tree(out_file, root, fields[1]);
    }
    std::cerr << "Done." << std::endl;
}

