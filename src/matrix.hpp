#ifndef MATRIX_HPP
#define MATRIX_HPP

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

class Matrix {
private:
    size_t default_row_number;                   // Default number of rows
    size_t column_number;                        // Number of columns
    std::vector<bool> matrix_1D;                 // we can try bitset instead 
    std::unordered_map<std::string, size_t> row_header;

public:
    // Constructor
    Matrix(size_t default_row_number = 1000000, size_t column_number = 2);

    // Getter for matrix
    const std::vector<bool>& get_matrix() const;

    // Setter for matrix
    void set_matrix(const std::vector<bool>& expanded_matrix);

    // Getter for row header
    const std::unordered_map<std::string, size_t>& get_row_header() const;

    // Getter for row count matrix
    size_t rowCount() const;

    // Getter for default row number
    size_t get_default_row_number() const;

    // expand matrix when is full (increase size by + default_row_number)
    void expandMatrix();

    // Setter for row header
    void set_row_header(const std::unordered_map<std::string, size_t>& row_header);

    // Add data at specified indices
    void add_data(size_t idx_snarl, size_t idx_geno);
};

#endif
