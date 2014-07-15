#include "Query.h"
#include "Build.h"
#include "BloomTree.h"
#include "BF.h"
#include "util.h"

#include <string>
#include <getopt.h>

/* TODO:
 * 2. Build a tree file from a collection of JF bloom filter
 *
 * (given leaves, build the 
 * 3. test the heap out for larger scale test
 * 4. make some constants be options
 */

// various commandline filenames
std::string command;
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
        << "Usage: bt [query|convert|build] ...\n"
        << "    \"query\" bloomtreefile queryfile\n"
        << "    \"convert\" jfbloomfilter outfile\n"
        << "    \"build\" filterlistfile outfile\n"
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
        if (optind >= argc-2) print_usage();
        bloom_tree_file = argv[optind+1];
        query_file = argv[optind+2];

    } else if (command == "convert") {
        if (optind >= argc-2) print_usage();
        jfbloom_file = argv[optind+1];
        out_file = argv[optind+2];

    } else if (command == "build") {
        if (optind >= argc-2) print_usage();
        query_file = argv[optind+1];
        bloom_tree_file = argv[optind+2];
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

        std::cerr << "Querying..." << std::endl;
        query_from_file(root, query_file, std::cout);

    } else if (command == "convert") {
        std::cerr << "Converting..." << std::endl;
        convert_jfbloom_to_rrr(jfbloom_file, out_file);

    } else if (command == "build") {
        std::cerr << "Building..." << std::endl;
        vector<std::string> leaves = read_filter_list(query_file);
        build_bt_from_jfbloom(leaves, out_file);
    }
    std::cerr << "Done." << std::endl;
}

