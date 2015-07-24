#ifndef QUERY_H
#define QUERY_H

#include <set>
#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "BloomTree.h"

extern float QUERY_THRESHOLD;

struct QueryInfo {
    QueryInfo(const std::string & q) : query(q), query_kmers(kmers_in_string(q)) {}
    QueryInfo(const std::string & q, const std::string & w){
	query = q;
	query_kmers = kmers_in_string(q);
	std::vector<std::string> fields;
	SplitString(w, ' ', fields);
	unsigned n = 0;
	for (const auto & w : fields){
		if (w!="") {
			// If string w has invalid letters after numbers, composite_string will return pointer (else null)
			// std::string::size_type not working here?
			std::size_t composite_string;
			//std::cerr << typeid(w).name() << std::endl;
			try{
				float value = std::stof(w, &composite_string);
	                        if (composite_string != w.size()) {
	                                std::cerr << "Invalid weight \'" << w << "\' at position " << n << std::endl;
        	                        exit(3);
                        	}
	                        weight.emplace_back(value); //Currently zero error handling here!

			}
			catch(...){
				std::cerr << "Invalid weight \'" << w << "\' at position " << n << std::endl;
				exit(3);
			}
		}
		n++;
    	}
    }
    ~QueryInfo() {}
   
    std::string query;
    std::set<jellyfish::mer_dna> query_kmers;
    std::vector<const BloomTree*> matching;
    std::vector<float> weight;
};

using QuerySet = std::list<QueryInfo*>;

void query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
void batch_weightedquery_from_file(BloomTree* root, const std::string & fn, const std::string & wf, std::ostream & o); 
void query_string(BloomTree* root, const std::string & q, std::vector<BloomTree*> & out);
void query(BloomTree* root, const std::set<jellyfish::mer_dna> & q, std::vector<BloomTree*> & out);
void check_bt(BloomTree* root);
void draw_bt(BloomTree* root, std::string outfile);
void compress_bt(BloomTree* root);

void leaf_query_from_file(BloomTree* root, const std::string & fn, std::ostream & o);
#endif
