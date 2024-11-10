#include "matrix.hpp"

// Constructor implementation
Matrix::Matrix(size_t default_row_number, size_t column_number)
    : default_row_number(default_row_number), column_number(column_number)
{
    unsigned long long int length_matrix = default_row_number * column_number; 
    matrix_1D.reserve(length_matrix);
    matrix_1D.resize(length_matrix, false);
}

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

void Matrix::expandMatrix() {
    unsigned long long int current_rows = matrix_1D.size() / column_number;  // Current number of rows
    unsigned long long int new_rows = current_rows + default_row_number;     // New total number of rows after expansion
    unsigned long long int new_matrix_size = new_rows * column_number;     // New total number of rows after expansion

    // Reserve space if needed, and resize the vector to the new size
    matrix_1D.reserve(new_matrix_size);          // Ensures capacity is sufficient
    matrix_1D.resize(new_matrix_size, false);    // Resizes the vector and fills new elements with 'false'
}

// Add data at specified indices
void Matrix::add_data(size_t idx_snarl, size_t idx_geno) {
    if (idx_snarl < matrix_1D.size() && idx_geno < column_number) {
        matrix_1D[idx_snarl * idx_geno] = true; // Set the value to true (1)
    } else {
        std::cerr << "Index out of bounds: (" << idx_snarl << ", " << idx_geno << ")" << std::endl;
    }
}
