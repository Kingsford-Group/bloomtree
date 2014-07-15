#ifndef BUILD_H
#define BUILD_H

#include <string>
#include <vector>

std::vector<std::string> read_filter_list(const std::string & inf);
void build_bloom_tree_filters(const std::vector<std::string> & fs, const std::string & matrix_file, const std::string & outf);
void convert_jfbloom_to_rrr(const std::string & jfbloom_file, const std::string & out_file);

#endif
