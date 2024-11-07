#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <chrono>
#include <boost/math/distributions/chi_squared.hpp>

// Function to calculate the Chi-square test statistic
double chiSquareStatistic(const std::vector<std::vector<int>>& observed) {
    int rows = observed.size();
    int cols = observed[0].size();

    if (rows < 2 || cols < 2) {
        throw std::invalid_argument("Input table must have at least 2 rows and 2 columns.");
    }

    // Calculate row and column sums
    std::vector<int> rowSums(rows, 0);
    std::vector<int> colSums(cols, 0);
    int totalSum = 0;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            rowSums[i] += observed[i][j];
            colSums[j] += observed[i][j];
            totalSum += observed[i][j];
        }
    }

    // Compute expected frequencies and Chi-square statistic
    double chiSquare = 0.0;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double expected = (double(rowSums[i]) * colSums[j]) / totalSum;
            chiSquare += std::pow(observed[i][j] - expected, 2) / expected;
        }
    }

    return chiSquare;
}

// Function to calculate the degrees of freedom for a 2D table
int calculateDegreesOfFreedom(int rows, int cols) {
    return (rows - 1) * (cols - 1);
}

// Function to perform the Chi-square test
std::string chi2Test(const std::vector<std::vector<int>>& observed) {
    // Ensure the table has at least 2 rows and 2 columns and all cells have non-zero counts
    int rows = observed.size();
    int cols = observed[0].size();

    if (rows >= 2 && cols >= 2) {
        // Check that all rows and columns sum to non-zero values
        bool valid = true;
        for (int i = 0; i < rows && valid; ++i) {
            valid &= std::accumulate(observed[i].begin(), observed[i].end(), 0) > 0;
        }
        for (int j = 0; j < cols && valid; ++j) {
            int colSum = 0;
            for (int i = 0; i < rows; ++i) {
                colSum += observed[i][j];
            }
            valid &= colSum > 0;
        }

        if (valid) {
            try {
                double chiSquare = chiSquareStatistic(observed);
                int degreesOfFreedom = calculateDegreesOfFreedom(rows, cols);

                // Use Boost to calculate the p-value from the Chi-square distribution
                boost::math::chi_squared dist(degreesOfFreedom);
                double pValue = 1 - cdf(dist, chiSquare); // 1 - CDF gives the p-value

                return std::to_string(pValue);  // Return the p-value as a string
            } catch (const std::exception& e) {
                return "Error: " + std::string(e.what());
            }
        } else {
            return "N/A";  // Invalid table data
        }
    } else {
        return "N/A";  // Invalid table dimensions
    }
}

// ------------------- Fisher -------------------

// Function to initialize the log factorials array
void initLogFacs(double* logFacs, int n) {
    logFacs[0] = 0.0;  // log(1) = 0, factorial of 0 is 1
    for (int i = 1; i <= n; ++i) {
        logFacs[i] = logFacs[i - 1] + log((double)i);  // log(n!) = log((n-1)!) + log(n)
    }
}

// Function to calculate the log probability of the hypergeometric distribution
double logHypergeometricProb(double* logFacs , int a, int b, int c, int d) {
    return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d]
    - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d];
}

double fastFishersExactTest(const std::vector<std::vector<int>>& table) {
    // Ensure the table is 2x2
    if (table.size() != 2 || table[0].size() != 2 || table[1].size() != 2) {
        throw std::invalid_argument("Input table must be 2x2.");
    }

    // Extract values from the table
    int a = table[0][0];
    int b = table[0][1];
    int c = table[1][0];
    int d = table[1][1];

    // Total sum of the table
    int n = a + b + c + d;

    // Dynamically allocate memory for logFacs array
    double* logFacs = new double[n + 1];
    if (!logFacs) {
        throw std::runtime_error("Memory allocation failed for logFacs.");
    }

    // Initialize log factorials
    initLogFacs(logFacs, n);

    // Compute log probability cutoff
    double logpCutoff = logHypergeometricProb(logFacs, a, b, c, d);
    double pFraction = 0.0;

    // Compute pFraction by iterating through possible values
    for (int x = 0; x <= n; ++x) {
        if (a + b - x >= 0 && a + c - x >= 0 && d - a + x >= 0) {
            double l = logHypergeometricProb(logFacs, x, a + b - x, a + c - x, d - a + x);
            if (l <= logpCutoff) pFraction += exp(l - logpCutoff);
        }
    }

    // Clean up memory
    delete[] logFacs;
    return exp(logpCutoff + log(pFraction));
}

int main() {
    std::vector<std::vector<int>> table = {{1982, 3018}, {2056, 2944}};

    // Measure time for Fisher's exact test
    clock_t start = clock();
    try {
        double p_value = fastFishersExactTest(table);
        std::cout << "Fisher's exact test p-value: " << p_value << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    clock_t end = clock();
    std::cout << "Time for fisherExactTest: " 
              << (double)(end - start) / CLOCKS_PER_SEC * 1e6  // Convert seconds to microseconds
              << " microseconds" << std::endl;

    // Measure time for Chi2 test
    auto start2 = std::chrono::high_resolution_clock::now();
    try {
        std::string p_value = chi2Test(table);
        std::cout << "Chi2 test p-value : " << p_value << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end2 - start2);
    std::cout << "Time for chi2Test: " << duration2.count() << " microseconds" << std::endl;

    return 0;
}
