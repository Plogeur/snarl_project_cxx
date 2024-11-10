#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <unordered_set>
#include <sstream>
#include <utility>
#include <iostream>
#include <chrono>
#include "matrix.hpp"

struct Variant {
    std::vector<std::string> atInfo;
    std::vector<std::vector<int>> genotypes;
    Variant(std::string& line, std::vector<std::string> sampleNames);
};

// SnarlParser class declaration
class SnarlParser {
public:
    Matrix matrix; // Matrix(size_t default_row_number = 1000000, size_t column_number = 2);
    const std::string& filename;
    std::vector<std::string> sampleNames;
    std::ifstream file;
    
    SnarlParser(const std::string& vcf_path);
    void pushMatrix(const std::string& decomposedSnarl, std::unordered_map<std::string, size_t>& rowHeaderDict, size_t indexColumn);
    void fill_matrix();
    std::vector<std::string> parseHeader();
    void binary_table(const std::unordered_map<std::string, std::vector<std::string> >& snarls,
                        const std::unordered_map<std::string, bool>& binary_groups,
                        const std::string& output = "output/binary_output.tsv");
    void quantitative_table(const std::unordered_map<std::string, std::vector<std::string> >& snarls,
                                const std::unordered_map<std::string, float>& quantitative,
                                const std::string& output = "output/quantitative_output.tsv");
};

// Retrieve the index of `key` if it exists in `ordered_map`. Otherwise, add it and return the new index.
size_t getOrAddIndex(std::unordered_map<std::string, size_t>& orderedMap, const std::string& key, size_t lengthOrderedMap);

// Function to decompose a string with snarl information
std::vector<std::string> decompose_string(const std::string& s);

// Function to split a string by a delimiter
std::vector<std::string> split(const std::string& str, char delimiter);

std::vector<int> extractGenotype(const std::string& genotypeStr);

std::vector<std::string> extractATField(const std::string& infoField);

std::vector<std::string> split(const std::string& s, char delimiter);

// Function to determine and extract an integer from the string
std::pair<size_t, size_t> determine_str(const std::string& s, size_t length_s, size_t i);

// Function to decompose a list of snarl strings
const std::vector<std::vector<std::string>> decompose_snarl(const std::vector<std::string>& lst);

std::vector<int> identify_correct_path(const std::vector<std::string>& decomposed_snarl, 
                                        const std::unordered_map<std::string, size_t>& row_headers_dict, 
                                        std::vector<int>& srr_save, const Matrix& matrix, const size_t num_cols);
