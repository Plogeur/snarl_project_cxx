#include "matrix.hpp"

// Constructor implementation
Matrix::Matrix(size_t default_row_number, size_t column_number)
    : default_row_number(default_row_number), column_number(column_number) {
    matrix_2D.resize(default_row_number, std::vector<bool>(column_number, false));
}

// Getter for matrix
const std::vector<std::vector<bool>>& Matrix::get_matrix() const {
    return matrix_2D;
}

// Getter for row header
const std::unordered_map<std::string, int>& Matrix::get_row_header() const {
    return row_header;
}

// Getter for default row number
size_t Matrix::get_default_row_number() const {
    return default_row_number;
}

// Getter for row number
size_t Matrix::rowCount() const {
    return matrix_2D.size();
}

// Setter for row header
void Matrix::set_row_header(const std::vector<std::string>& row_header) {
    this->row_header = row_header;
}

// Setter for matrix
void Matrix::set_matrix(const std::vector<std::vector<bool>>& expanded_matrix) {
    this->matrix_2D = expanded_matrix;
}

// Method to expand the matrix by adding rows
void Matrix::expandMatrix() {
    int current_rows = matrix_2D.size();  // Number of rows
    int current_cols = current_rows > 0 ? matrix_2D[0].size() : 0;  // Number of columns (if rows are non-zero)
    int new_rows = current_rows + default_row_number;

    // Create a new matrix with the expanded size
    std::vector<std::vector<bool>> expanded_matrix(new_rows, std::vector<bool>(current_cols, false));

    // Copy the current matrix data to the expanded matrix
    for (int i = 0; i < current_rows; ++i) {
        for (int j = 0; j < current_cols; ++j) {
            expanded_matrix[i][j] = matrix_2D[i][j];
        }
    }

    // Update the matrix data with the expanded matrix
    set_matrix(expanded_matrix);
}

// Add data at specified indices
void Matrix::add_data(size_t idx_snarl, size_t idx_geno) {
    if (idx_snarl < matrix_2D.size() && idx_geno < column_number) {
        matrix_2D[idx_snarl][idx_geno] = true; // Set the value to true (1)
    } else {
        std::cerr << "Index out of bounds: (" << idx_snarl << ", " << idx_geno << ")" << std::endl;
    }
}
