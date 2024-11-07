#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <chrono>  // Include the chrono library for timing

// Helper function for summing over axes (in this case, the entire vector)
std::vector<double> apply_over_axes_sum(const std::vector<std::vector<double>>& a, const std::vector<int>& axes_to_sum) {
    std::vector<double> result(a.size(), 0.0);

    // Sum along the specified axes
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < a[i].size(); ++j) {
            result[i] += a[i][j];
        }
    }

    return result;
}

// Marginal sums function (equivalent to Python's margins)
std::vector<std::vector<double>> margins(const std::vector<std::vector<double>>& a) {
    std::vector<std::vector<double>> margsums;
    size_t rows = a.size();
    size_t cols = a[0].size();

    for (size_t k = 0; k < cols; ++k) {
        std::vector<int> remaining_axes;
        for (size_t i = 0; i < cols; ++i) {
            if (i != k) remaining_axes.push_back(i);
        }

        std::vector<double> marg = apply_over_axes_sum(a, remaining_axes);
        margsums.push_back(marg);
    }

    return margsums;
}

// Compute the expected frequencies
std::vector<std::vector<double>> expected_freq(const std::vector<std::vector<double>>& observed) {
    size_t rows = observed.size();
    size_t cols = observed[0].size();

    // Get marginal sums
    std::vector<std::vector<double>> margsums = margins(observed);

    // Perform the computation for expected frequencies
    double total_sum = 0;
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            total_sum += observed[i][j];
        }
    }

    std::vector<std::vector<double>> expected(rows, std::vector<double>(cols, 0.0));
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            expected[i][j] = (margsums[i][j] * margsums[j][i]) / std::pow(total_sum, 2);
        }
    }
    return expected;
}

// Power divergence (a simplified placeholder)
double power_divergence(const std::vector<std::vector<double>>& observed, const std::vector<std::vector<double>>& expected) {
    double sum = 0.0;

    for (size_t i = 0; i < observed.size(); ++i) {
        for (size_t j = 0; j < observed[0].size(); ++j) {
            double diff = observed[i][j] - expected[i][j];
            sum += diff * diff; // Simplified squared difference
        }
    }

    return sum;
}

// Chi-squared test calculation
std::tuple<double, double, int, std::vector<std::vector<double>>> chi2_contingency(std::vector<std::vector<double>>& observed, bool correction=true) {
    // Check for negative values in observed data
    for (size_t i = 0; i < observed.size(); ++i) {
        for (size_t j = 0; j < observed[0].size(); ++j) {
            if (observed[i][j] < 0) {
                throw std::invalid_argument("All values in `observed` must be nonnegative.");
            }
        }
    }

    if (observed.empty()) {
        throw std::invalid_argument("No data; `observed` has size 0.");
    }

    // Compute expected frequencies
    std::vector<std::vector<double>> expected = expected_freq(observed);

    // Check if any expected value is zero
    for (size_t i = 0; i < expected.size(); ++i) {
        for (size_t j = 0; j < expected[0].size(); ++j) {
            if (expected[i][j] == 0) {
                throw std::invalid_argument("The internally computed table of expected frequencies has a zero element.");
            }
        }
    }

    // Calculate degrees of freedom
    int dof = expected.size() * expected[0].size() - observed.size() - observed[0].size() + 1;

    // In case of a degenerate case
    double chi2 = 0.0;
    double p = 1.0;
    if (dof == 0) {
        chi2 = 0.0;
        p = 1.0;
    } else {
        // Apply Yates' correction if needed
        if (correction) {
            for (size_t i = 0; i < observed.size(); ++i) {
                for (size_t j = 0; j < observed[0].size(); ++j) {
                    double diff = expected[i][j] - observed[i][j];
                    double direction = (diff > 0) ? 1.0 : -1.0;
                    double magnitude = std::min(0.5, std::abs(diff));
                    observed[i][j] += magnitude * direction;
                }
            }
        }

        // Calculate the chi-squared statistic
        chi2 = power_divergence(observed, expected);
        p = 1.0 - std::exp(-chi2 / 2);  // Placeholder for p-value calculation
    }

    return std::make_tuple(chi2, p, dof, expected);
}

int main() {
    // Example usage
    std::vector<std::vector<double>> table = {{1982, 3018}, {2056, 2944}};

    // Start measuring time
    auto start = std::chrono::high_resolution_clock::now();

    try {
        auto [chi2_stat, p_value, dof, expected_freq] = chi2_contingency(table);

        // End measuring time
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;

        std::cout << "Chi2 Statistic: " << chi2_stat << std::endl;
        std::cout << "p-value: " << p_value << std::endl;
        std::cout << "Degrees of Freedom: " << dof << std::endl;

        // Output the time taken
        std::cout << "Time taken for chi-squared calculation: " << duration.count() << " seconds." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
