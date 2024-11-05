#ifndef SNARL_PARSER_HPP
#define SNARL_PARSER_HPP

#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <utility> // For std::pair
#include <iostream> // For std::cout, std::cerr

// Function to split a string by a delimiter
std::vector<std::string> split(const std::string& str, char delimiter);

// Function to determine and extract an integer from the string
std::pair<int, int> determine_str(const std::string& s, int length_s, int i);

// Function to decompose a string with snarl information
std::vector<std::string> decompose_string(const std::string& s);

// Function to decompose a list of snarl strings
std::vector<std::vector<std::string>> decompose_snarl(const std::vector<std::string>& lst);

// Main function that parses the VCF file and fills the matrix
void fill_matrix(const std::string& vcf_path);

#endif
