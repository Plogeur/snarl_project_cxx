#include "binary_analysis.hpp"
#include "snarl_parser.hpp"

// ------------------------ Chi2 exact test ------------------------

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

                return std::to_string(pValue);
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

// ------------------------ Fisher exact test ------------------------

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

std::vector<std::vector<int>> create_binary_table(
    const std::unordered_map<std::string, bool>& groups, 
    const std::vector<std::string>& list_path_snarl) 
{
    std::map<std::string, int> row_headers_dict = matrix.get_row_header();  // Assuming matrix is accessible
    std::vector<std::string> list_samples = this->list_samples;             // Assuming list_samples is accessible
    size_t length_column_headers = list_path_snarl.size();

    // Initialize g0 and g1 with zeros, corresponding to the length of column_headers
    std::vector<int> g0(length_column_headers, 0);
    std::vector<int> g1(length_column_headers, 0);

    // Iterate over each path_snarl in column_headers
    for (size_t idx_g = 0; idx_g < list_path_snarl.size(); ++idx_g) {
        const std::string& path_snarl = list_path_snarl[idx_g];
        std::vector<int> idx_srr_save(list_samples.size());
        std::iota(idx_srr_save.begin(), idx_srr_save.end(), 0);  // Fill with [0, 1, ..., len(list_samples) - 1]

        std::vector<std::string> decomposed_snarl = decompose_string(path_snarl);
        idx_srr_save = identify_correct_path(decomposed_snarl, row_headers_dict, idx_srr_save);

        // Count occurrences in g0 and g1 based on the updated idx_srr_save
        for (int idx : idx_srr_save) {
            std::string srr = list_samples[idx / 2];  // Convert index to the appropriate sample name

            // Check if sample belongs to a group and increment the respective count
            auto it = groups.find(srr);
            if (it != groups.end()) {
                if (it->second) {
                    // If true, consider as part of Group 1
                    g1[idx_g] += 1;
                } else {
                    // If false, consider as part of Group 0
                    g0[idx_g] += 1;
                }
            } else {
                throw std::runtime_error("Sample " + srr + " not found in groups.");
            }
        }
    }

    // Return the populated binary table
    return {g0, g1};
}
