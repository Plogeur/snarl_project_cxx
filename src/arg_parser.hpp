#ifndef ARG_PARSERS_HPP
#define ARG_PARSERS_HPP

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <limits>
#include <filesystem>
#include <stdexcept>

// Parses the group file and fills the group_0 and group_1 maps with sample data.
std::unordered_map<std::string, bool> parse_group_file(const std::string& group_file);

// Parses the phenotype file and returns a map with IID as keys and PHENO as float values.
std::unordered_map<std::string, float> parse_pheno_file(const std::string& file_path);

// Parses the snarl path file and returns a map with snarl as keys and paths as a list of strings.
std::unordered_map<std::string, std::vector<std::string>> parse_snarl_path_file(const std::string& path_file);

std::string check_format_vcf_file(const std::string& file_path);
std::string check_format_group_snarl(const std::string& file_path);
std::string check_format_pheno_q(const std::string& file_path);
std::string check_format_pheno_b(const std::string& file_path);

#endif
