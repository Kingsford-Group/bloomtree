#include "Query.h"
#include "Build.h"
#include "BloomTree.h"
#include "BF.h"
#include "util.h"

#include <string>
#include <cstdlib>
#include <getopt.h>

/* TODO:
 */

// various commandline filenames
std::string command;
std::string bloom_tree_file;
std::string query_file;
std::string out_file;
std::string jfbloom_file;
std::string bvfile1, bvfile2;

std::string bloom_storage;

unsigned parallel_level = 3; // no parallelism by default

const char * OPTIONS = "t:p:f:";

static struct option LONG_OPTIONS[] = {
    {"max-filters", required_argument, 0, 'f'},
    {"threads", required_argument, 0, 'p'},
    {"query-threshold", required_argument, 0, 't'},
    {0,0,0,0}
};

void print_usage() {
    std::cerr 
        << "Usage: bt [query|convert|build] ...\n"
        << "    \"query\" [--max-filters 32] [-t 0.8] bloomtreefile queryfile outfile\n"
        << "    \"convert\" jfbloomfilter outfile\n"
        << "    \"build\" filterlistfile outfile\n"
        << "    \"check\" bloomtreefile\n"
        << "    \"sim\" bvfile1 bvfile2\n"
        << "    \"draw\" bloomtreefile out.dot\n"
        << std::endl;
    exit(3);
}

int process_options(int argc, char* argv[]) {
    int a;
    while ((a=getopt_long(argc, argv, OPTIONS, LONG_OPTIONS, 0)) != -1) {
        switch(a) {
            case 't':
                QUERY_THRESHOLD = atof(optarg);
                break;
            case 'p':
                parallel_level = unsigned(atoi(optarg));
                break;
            case 'f':
                BF_INMEM_LIMIT = unsigned(atoi(optarg));
                break;
            default:
                std::cerr << "Unknown option." << std::endl;
                print_usage();
        }
    }

    if (optind >= argc) print_usage();
    command = argv[optind];
    if (command == "query") {
        if (optind >= argc-3) print_usage();
        bloom_tree_file = argv[optind+1];
        query_file = argv[optind+2];
	out_file = argv[optind+3];
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
    } else if (command == "convert") {
        if (optind >= argc-2) print_usage();
        jfbloom_file = argv[optind+1];
        out_file = argv[optind+2];
    } else if (command == "build") {
        if (optind >= argc-3) print_usage();
        query_file = argv[optind+1];
        out_file = argv[optind+2];
        bloom_storage = argv[optind+3];
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
        batch_query_from_file(root, query_file, out);

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
        std::tuple<uint64_t, uint64_t> sim = bf1->b_similarity(bf2);
	std::cerr << "Done " << std::endl;
	std::cout << bf1->size() << " " << std::get<0>(sim) << " " << std::get<1>(sim) << std::endl;

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

    } else if (command == "build") {
        std::cerr << "Building..." << std::endl;
        std::vector<std::string> leaves = read_filter_list(query_file);
        //build_bt_from_jfbloom(leaves, out_file, parallel_level);
        dynamic_build(leaves, out_file, bloom_storage);
    }
    std::cerr << "Done." << std::endl;
}

