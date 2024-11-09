#include "matrix.hpp"

// Constructor implementation
Matrix::Matrix(size_t default_row_number, size_t column_number)
    : default_row_number(default_row_number), column_number(column_number),
    matrix_1D(default_row_number * column_number, false)
{}

// Getter for matrix
const std::vector<bool>& Matrix::get_matrix() const {
    return matrix_1D;
}

// Getter for row header
const std::unordered_map<std::string, size_t>& Matrix::get_row_header() const {
    return row_header;
}

// Getter for default row number
size_t Matrix::get_default_row_number() const {
    return default_row_number;
}

// Getter for row number
size_t Matrix::rowCount() const {
    return matrix_1D.size();
}

// Setter for row header
void Matrix::set_row_header(const std::unordered_map<std::string, size_t>& row_header) {
    this->row_header = row_header;
}

// Setter for matrix
void Matrix::set_matrix(const std::vector<bool>& expanded_matrix) {
    this->matrix_1D = expanded_matrix;
}

// Method to expand the matrix by adding rows
void Matrix::expandMatrix() {
    int current_rows = matrix_1D.size()/column_number;  // Number of rows
    int new_rows = current_rows + default_row_number;

    // Create a new matrix with the expanded size
    std::vector<bool> expanded_matrix(new_rows * column_number, false);
    expanded_matrix = matrix_1D;

    // Update the matrix data with the expanded matrix
    set_matrix(expanded_matrix);
}

// Add data at specified indices
void Matrix::add_data(size_t idx_snarl, size_t idx_geno) {
    if (idx_snarl < matrix_1D.size() && idx_geno < column_number) {
        matrix_1D[idx_snarl * idx_geno] = true; // Set the value to true (1)
    } else {
        std::cerr << "Index out of bounds: (" << idx_snarl << ", " << idx_geno << ")" << std::endl;
    }
}
