#include "Build.h"
#include "BloomTree.h"
#include "util.h"
#include "Query.h"
#include <cmath>
#include <sstream>
#include <cstring>
#include <future>
#include <algorithm>
#include "gzstream.h"

#include <sys/mman.h>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/file_header.hpp>
#include <sdsl/bit_vectors.hpp>

// non-zero if bit is set
inline char bit(char * buf, unsigned long bit) {
    char byte = buf[bit / 8];
    char bit_mask = 1 << (bit % 8); 
    return byte & bit_mask;
}


sdsl::bit_vector* read_bit_vector_from_jf(const std::string & jfbloom_file) {
    // Load the JF Bloom Filter
    std::cerr << "Reading: " << jfbloom_file << std::endl;
    igzstream in(jfbloom_file.c_str(), std::ios::in|std::ios::binary);
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
    sdsl::bit_vector* b = new sdsl::bit_vector(length, 0);
    std::memcpy(b->data(), buf, num_bytes);
    //std::memcpy((void*) b->data(), buf, num_bytes);

    /*for (unsigned long i = 0; i < length; i++) {
        if (bit(buf, i)) {
            (*b)[i] = 1;
        }
    }*/
    delete buf;
    return b;
}


void convert_jfbloom_to_rrr(const std::string & jfbloom_file, const std::string & out_file) {
    sdsl::bit_vector* b = read_bit_vector_from_jf(jfbloom_file);
    
    // covert that raw bit vector into an rrr compressed vector & save it.
    sdsl::rrr_vector<255> rrr(*b); 
    std::cerr << "Compressed RRR vector is " << sdsl::size_in_mega_bytes(rrr) << std::endl; 
    sdsl::store_to_file(rrr, out_file);

    delete b;
}


// read a list of bloom filters, 1 per line from the given file
std::vector<std::string> read_filter_list(const std::string & inf) {
    std::ifstream in(inf);
    std::vector<std::string> v;
    std::string line;
    while(getline(in, line)) {
        v.emplace_back(Trim(line));
    }
    std::cerr << "Building tree on " << v.size() << " bloom filters." << std::endl;
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


// compute the number of nodes needed for a full binary tree of n nodes
unsigned number_nodes_in_complete_tree(unsigned n) {
    unsigned lp2 = int(log2(n));  // largest power of 2 <= n
    unsigned sz_last_complete_row = pow(2, lp2);
    return (2 * sz_last_complete_row - 1) + 2*(n-sz_last_complete_row);
}


/*
// a very simple implementation of union
sdsl::bit_vector* union_bv(const sdsl::bit_vector& b1, const sdsl::bit_vector& b2) {
    sdsl::bit_vector* out = new sdsl::bit_vector(b1.size(), 0);
    for (std::size_t i = 0; i < b1.size(); i++) {
        (*out)[i] = b1[i] | b2[i];
    }
    return out;
}
*/



// do a post-order traversal over the array-based "tree"
sdsl::bit_vector* build_filters(
    const std::vector<std::string> & leaves,
    std::vector<BloomTree*> & tree,
    std::vector<sdsl::bit_vector*> & raw,
    const HashPair & hashes, 
    int nh,
    std::size_t pos
) {
    // build the left bitvector (and saved rrr vector)
    sdsl::bit_vector* bl = nullptr;
    auto left = complete_tree_child(pos, 0);
    if (left < tree.size()) {
        if (raw[left] == nullptr) {
            bl = build_filters(leaves, tree, raw, hashes, nh, left);
        } else {
            bl = raw[left];
        }
    }

    // build the right bitvector (and saved rrr vector)
    sdsl::bit_vector* br = nullptr;
    auto right = complete_tree_child(pos, 1);
    if (right < tree.size()) {
        if (raw[right] == nullptr) {
            br = build_filters(leaves, tree, raw, hashes, nh, right);
        } else {
            br = raw[right];
        }
    }

    sdsl::bit_vector* u = nullptr;
    std::string union_name;

    // if interior node:
    if (br != nullptr && bl != nullptr) {
        std::cerr << "Creating union: " << pos << std::endl;
        
        // union the two filters & discard the children
        u = union_bv_fast(*bl, *br);
        raw[pos] = u;
        raw[left] = nullptr;
        raw[right] = nullptr;
        delete bl;
        delete br;
        
        // create the new name
        std::ostringstream oss;
        oss << "union" << pos << ".rrr";
        union_name = oss.str();

        std::cerr << "Unioning: " << tree[left]->name() << " with " <<
            tree[right]->name() << " to " << union_name << std::endl;

        // save the BT node
        tree[pos] = new BloomTree(union_name, hashes, nh);
        tree[pos]->set_child(0, tree[left]);
        tree[pos]->set_child(1, tree[right]);

    } else if (br == nullptr && br == nullptr) {
        // we're a leaf, so what we should do is (1) read the JF bloom filter,
        // (2) create a bv, (3) store a compressed rrr vector
        union_name = leaves[tree.size() - pos - 1];
        u = read_bit_vector_from_jf(union_name);
        raw[pos] = u;

        union_name = test_basename(union_name, std::string(".gz")) + ".rrr";
        tree[pos] = new BloomTree(union_name, hashes, nh);

    } else {
        DIE("Should not happen.");
    }
    
    // compress and write it out if it doesn't already exist
    std::ifstream check(union_name.c_str());
    if (!check) {
        // convert and save compressed version
        std::cerr << "Compressing to " << union_name << std::endl;
        sdsl::rrr_vector<255> rrr(*u);
        sdsl::store_to_file(rrr, union_name);
        std::cerr << "Compressed RRR vector is " << sdsl::size_in_mega_bytes(rrr) << std::endl; 
    } else {
        check.close();
        std::cerr << "Skipping compression because " << union_name << " already exists." << std::endl;
    }

    return u;
}


sdsl::bit_vector* build_filters_parallel(
    const std::vector<std::string> & leaves, 
    std::vector<BloomTree*> & tree,
    const HashPair & hashes, 
    int nh, 
    unsigned level
) {
    std::vector<sdsl::bit_vector*> raw(tree.size(), nullptr);

    /*
    level 1 = [1,2]
    level 2 = [3,4,5,6]
    level 3 = [7,8,9,10,11,12,13, 14]
    level 4 = [15,16,17,18, 19,20,21,22, 23,24,25,26, 27,28,29,30]
    level i starts at 2^i-1, and continues for 2^i
    */

    unsigned biggest_complete_level = unsigned(log2(leaves.size()));
    level = std::min(level, biggest_complete_level);

    // the roots of the subtrees that we are going to run in parallel
    std::size_t level_length = std::size_t(pow(2, level));
    std::size_t level_start = level_length - 1;

    std::cerr << "Parallel level = " << level << std::endl;
    std::cerr << "Computing " << level_length << " roots at level " << level <<
        " in parallel." << std::endl; 

    // for each of the roots, start a build_filters() process.
    std::vector<std::future<sdsl::bit_vector*> > F;
    for (std::size_t p = level_start; p < level_start + level_length; ++p) {
        F.emplace_back(
            std::async(
                std::launch::async,
                [=, &leaves, &tree, &raw] () {
                    return build_filters(leaves, tree, raw, hashes, nh, p);
                }
            )
        );
    }

    std::cerr << "Waiting for subtrees to finish." << std::endl;
    // join on all the threads
    for (auto & f : F) {
        f.wait();
    }

    // finish off the top of the tree
    std::cerr << "Finishing top of tree in single thread mode." << std::endl;
    return build_filters(leaves, tree, raw, hashes, nh, 0);
}


void build_bt_from_jfbloom(
    const std::vector<std::string> & leaves, 
    const std::string & outf,
    unsigned parallel_level
) {
    // create the hashes
    int nh = 0;
    HashPair* hashes = get_hash_function(leaves[0], nh); 

    // calculate the number of nodes in the tree
    unsigned nb_nodes = number_nodes_in_complete_tree(leaves.size());
    std::cerr << "Tree will have " << nb_nodes << " nodes" << std::endl;

    // v holds the nodes of the semi-complete tree
    std::vector<BloomTree*> v(nb_nodes, nullptr);
    sdsl::bit_vector *u = build_filters_parallel(leaves, v, *hashes, nh, parallel_level);
    std::cerr << "Built the whole tree." << std::endl;
    write_bloom_tree(outf, v[0], leaves[0]);

    // save uncompressed root for later merging
    std::cerr << "Storing the root bitvector" << std::endl;
    sdsl::store_to_file(*u, outf + ".root");
    
    // free the memory for the bloom nodes
    delete u;
    for (auto & p : v) {
        delete p;
    }
}


// walk down T, finding the best path; insert N (which could be a subtree) at the leaf
// we come to, and union all the parents
BloomTree* insert_bloom_tree(BloomTree* T, BloomTree* N, int type) {

    std::cerr << "Inserting leaf " << N->name() << " ... " << std::endl;
    // save the root to return
    BloomTree* root = T;

    // handle the case of inserting into an empty tree
    if (T == nullptr) {
        N->set_parent(nullptr);
        return N;
    }

    // until we fall off the tree (should insert before then)
    int depth = 0;
    int best_child = -1;
    BloomTree* parent = nullptr;

    while (T != nullptr) {
        std::cerr << "At node: " << T->name() << std::endl;
        T->increment_usage();
        if (T->num_children() == 0) {
            // this is the tricky case: T is currently a leaf, which means it
            // represents an SRA file, and so it has to stay a leaf. So what we
            // must do is replace T by a new union fiilter T -->
            // NewNode{child0=T, child1=N}
            std::ostringstream oss;
            oss << nosuffix(N->name(), std::string(".bf.bv")) << "_union.bf.bv";
            std::cerr << "Splitting leaf into " << oss.str() 
                << " at depth " << depth << std::endl;

            BloomTree* NewNode = T->union_bloom_filters(oss.str(), N);
            std::cerr << "   1:" << NewNode->child(0)->name() << std::endl;
            std::cerr << "   2:" << NewNode->child(1)->name() << std::endl;

            if (parent == nullptr) {
                return NewNode;
            } else {
                assert(best_child != -1);
                assert(parent != nullptr);
                parent->set_child(best_child, NewNode);
                return root;
            }
        } else if (T->num_children() == 1) {
            // union the new filter with this node
            T->union_into(N);
            
            // insert into first empty child
            for (int i =0; i < 2; i++) {
                if (T->child(i) == nullptr) {
                    std::cerr << "Adding as " << ((i==0)?"left":"right") 
                        << " child." << std::endl;
                    T->set_child(i, N);
                    return root;
                }
            }
            DIE("Something is wrong!");
        } else {
            // find the most similar child and move to it
            uint64_t best_sim = 0;
            best_child = -1;
            for (int i = 0; i < 2; i++) {
                T->child(i)->increment_usage();
                uint64_t sim = T->child(i)->similarity(N,type);
                std::cerr << "Child " << i << " sim =" << sim << std::endl;
                if (sim >= best_sim) {
                    best_sim = sim;
                    best_child = i;
                }
            }
            
            // union the new filter with this node
            T->union_into(N);

            // move the current ptr to the most similar child
            std::cerr << "Moving to " << ((best_child==0)?"left":"right") 
                << " child: " << best_child << " " << best_sim << std::endl;
            parent = T;
            T = T->child(best_child);
        }
        depth++;
    }
    assert(T != 0);

    return root;
}

// delete a tree
void delete_bloom_tree(BloomTree* T) {
    if (T->child(0) != nullptr) {
        delete_bloom_tree(T->child(0));
    }
    if (T->child(1) != nullptr) {
        delete_bloom_tree(T->child(1));
    }
    delete T;
}

// build the tree by repeated insertion
void dynamic_build(
    const std::string & hashes_file,
    const std::vector<std::string> & leaves, 
    const std::string & outf,
    //const std::string & bloom_storage,
    int type
) {
    // create the hashes
    int nh = 0;
    HashPair* hashes = get_hash_function(hashes_file, nh); 

    BloomTree* root = nullptr;
    
    //int count = 0;
    // for every leaf
    std::cerr << "Inserting leaves into tree..." << std::endl;
    for (const auto & leaf : leaves) {
        /*
        // read JF filter and save it as a bv
        sdsl::bit_vector* f = read_bit_vector_from_jf(leaf);
        std::string filter_name = test_basename(leaf, std::string(".gz")) + ".bv";
        std::string store_name = bloom_storage;
        store_name.append(filter_name);
        sdsl::store_to_file(*f, store_name);
        */

        // create the node that points to the filter we just saved
        BloomTree* N = new BloomTree(leaf, *hashes, nh);
        
        //  insert this new leaf
        root = insert_bloom_tree(root, N, type);

	/*
        count++;
        if (count % 100 == 0){
            std::string temp_out = outf;
            temp_out.append("_");
            temp_out.append(std::to_string(count));
            draw_bt(root, temp_out);
        }
	*/
    }
    
    // write the tree file
    std::cerr << "Built the whole tree." << std::endl;
    write_bloom_tree(outf, root, hashes_file);
    
    // delete the tree (which saves it)
    std::cerr << "Freeing tree (and saving dirty filters)" << std::endl;
    delete_bloom_tree(root);
}


/*** NOT USED ***/
/*
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
*/
