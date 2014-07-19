#ifndef BUILD_H
#define BUILD_H

#include <string>
#include <vector>

std::vector<std::string> read_filter_list(const std::string & inf);
void convert_jfbloom_to_rrr(const std::string & jfbloom_file, const std::string & out_file);
void build_bt_from_jfbloom(const std::vector<std::string> & leaves, const std::string & outf, unsigned parallel_level);

#endif
