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
    std::vector<std::vector<bool> > matrix_2D;      // 2D vector to represent the matrix
    std::unordered_map<std::string, int> row_header;

public:
    // Constructor
    Matrix(size_t default_row_number = 1000000, size_t column_number = 2);

    // Getter for matrix
    const std::vector<std::vector<bool> >& get_matrix() const;

    // Setter for matrix
    void set_matrix(const std::vector<std::vector<bool> >& expanded_matrix);

    // Getter for row header
    const std::unordered_map<std::string, int>& get_row_header() const;

    size_t rowCount() const;

    // Getter for default row number
    size_t get_default_row_number() const;

    void expandMatrix();

    // Setter for row header
    void set_row_header(const std::unordered_map<std::string, int>& row_header);

    // Add data at specified indices
    void add_data(size_t idx_snarl, size_t idx_geno);
};

#endif // MATRIX_HPP
