#include "Matrix.hpp"

// Constructor implementation
Matrix::Matrix(size_t default_row_number, size_t column_number)
    : default_row_number(default_row_number), column_number(column_number) {
    // Initialize the matrix with the default size filled with false (0)
    matrix.resize(default_row_number, std::vector<bool>(column_number, false));
}

// Getter for matrix
const std::vector<std::vector<bool>>& Matrix::get_matrix() const {
    return matrix;
}

// Setter for matrix
void Matrix::set_matrix(const std::vector<std::vector<bool>>& expanded_matrix) {
    matrix = expanded_matrix;
}

// Getter for row header
const std::vector<std::string>& Matrix::get_row_header() const {
    return row_header;
}

// Getter for default row number
size_t Matrix::get_default_row_number() const {
    return default_row_number;
}

// Getter for row number
size_t Matrix::rowCount() const {
    return matrix.size();
}

// Setter for row header
void Matrix::set_row_header(const std::vector<std::string>& row_header) {
    this->row_header = row_header;
}

// Add data at specified indices
void Matrix::add_data(size_t idx_snarl, size_t idx_geno) {
    if (idx_snarl < matrix.size() && idx_geno < column_number) {
        matrix[idx_snarl][idx_geno] = true; // Set the value to true (1)
    } else {
        std::cerr << "Index out of bounds: (" << idx_snarl << ", " << idx_geno << ")" << std::endl;
    }
}
