#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <unordered_set>
#include <sstream>
#include <utility> // For std::pair
#include <iostream> // For std::cout, std::cerr
#include "matrix.hpp"

// SnarlParser class declaration
class SnarlParser {
public:
    SnarlParser(const std::string& vcf_path);

    // Function to split a string by a delimiter
    std::vector<std::string> split(const std::string& str, char delimiter);

    // Function to determine and extract an integer from the string
    std::pair<int, int> determine_str(const std::string& s, int length_s, int i);

    // Function to decompose a string with snarl information
    std::vector<std::string> decompose_string(const std::string& s);

    // Function to decompose a list of snarl strings
    std::vector<std::vector<std::string> > decompose_snarl(const std::vector<std::string>& lst);
    
    // Retrieve the index of `key` if it exists in `ordered_map`. Otherwise, add it and return the new index.
    size_t getOrAddIndex(std::unordered_map<std::string, int>& orderedMap, const std::string& key, int lengthOrderedMap);

    // Add True to the matrix if snarl is found
    void pushMatrix(const std::string& decomposedSnarl, std::unordered_map<std::string, int>& rowHeaderDict, size_t indexColumn);

    // Main function that parses the VCF file and fills the matrix
    void fill_matrix();

    void binary_table(const std::unordered_map<std::string, std::vector<std::string> >& snarls,
                        const std::unordered_map<std::string, bool>& binary_groups,
                        const std::string& output = "output/binary_output.tsv");

    void quantitative_table(const std::unordered_map<std::string, std::vector<std::string> >& snarls,
                                const std::unordered_map<std::string, float>& quantitative,
                                const std::string& output = "output/quantitative_output.tsv");

private:
    std::vector<std::string> list_samples;
    Matrix matrix; // Matrix(size_t default_row_number = 1000000, size_t column_number = 2);
    std::string vcf_path;
};
