#include "quantitative_analysis.hpp"
#include "snarl_parser.hpp"

// Function to calculate p-values based on the t-statistic
double calculatePValue(double t_stat, int degrees_of_freedom) {
    // Approximation using the complementary error function
    double p_value = 2 * (1 - std::erf(std::fabs(t_stat) / std::sqrt(2)));
    return p_value;
}

void linearRegression(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Eigen::VectorXd& beta, Eigen::VectorXd& se, Eigen::VectorXd& p_values) {
    int n = X.rows();
    int p = X.cols();

    // Calculate beta using the OLS formula: beta = (X'X)^-1 X'y
    Eigen::MatrixXd XtX_inv = (X.transpose() * X).inverse();
    beta = XtX_inv * X.transpose() * y;

    // Calculate residuals and standard error
    Eigen::VectorXd y_pred = X * beta;
    Eigen::VectorXd residuals = y - y_pred;

    // Residual sum of squares (RSS)
    double rss = residuals.squaredNorm();
    double sigma2 = rss / (n - p);

    // Standard error of beta coefficients: SE(beta) = sqrt(sigma^2 * diag((X'X)^-1))
    se = (sigma2 * XtX_inv.diagonal()).array().sqrt();

    // Calculate t-statistics for beta coefficients
    Eigen::VectorXd t_stats = beta.array() / se.array();

    // Calculate p-values for each coefficient
    p_values.resize(p);
    for (int i = 0; i < p; ++i) {
        p_values[i] = calculatePValue(t_stats[i], n - p);
    }
}

std::vector<std::vector<int>> create_quantitative_table(const std::vector<std::string>& column_headers) {
    std::map<std::string, int> row_headers_dict = matrix.get_row_header();
    int length_sample = list_samples.size();

    // Initialize a zero matrix for genotypes with shape (length_sample, len(column_headers))
    std::vector<std::vector<int>> genotypes(length_sample, std::vector<int>(column_headers.size(), 0));

    // Iterate over each path_snarl and fill in the matrix
    for (size_t col_idx = 0; col_idx < column_headers.size(); ++col_idx) {
        const std::string& path_snarl = column_headers[col_idx];
        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);

        // Identify correct paths
        std::vector<int> idx_srr_save = identify_correct_path(decomposed_snarl, row_headers_dict, std::vector<int>(length_sample));

        for (int idx : idx_srr_save) {
            int srr_idx = idx / 2;  // Convert index to the appropriate sample index
            genotypes[srr_idx][col_idx] += 1;
        }
    }

    // Return the populated genotypes matrix
    return genotypes;
}


