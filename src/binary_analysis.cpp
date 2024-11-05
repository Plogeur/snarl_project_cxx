#include "binary_analysis.hpp"
#include <iostream>
#include <cmath>   // For sqrt, exp, and log
#include <boost/math/distributions/chi_squared.hpp>

// Function to calculate the chi-squared statistic and p-value
double chiSquaredTest(const Eigen::MatrixXd& observed, double& p_value) {
    Eigen::MatrixXd expected = observed.colwise().sum() * observed.rowwise().sum().transpose() / observed.sum();
    double chi_squared = ((observed - expected).array().square() / expected.array()).sum();

    // Calculate degrees of freedom
    int df = (observed.rows() - 1) * (observed.cols() - 1);

    // Calculate p-value using the chi-squared distribution
    boost::math::chi_squared chi_squared_dist(df);
    p_value = 1 - boost::math::cdf(chi_squared_dist, chi_squared);

    return chi_squared;
}

// Function to calculate Fisher's Exact Test p-value
// This is a simple implementation for a 2x2 contingency table
double fisherExactTest(const Eigen::MatrixXd& contingencyTable) {
    if (contingencyTable.rows() != 2 || contingencyTable.cols() != 2) {
        throw std::invalid_argument("Fisher's exact test only supports 2x2 contingency tables.");
    }

    double a = contingencyTable(0, 0);
    double b = contingencyTable(0, 1);
    double c = contingencyTable(1, 0);
    double d = contingencyTable(1, 1);
    double n = a + b + c + d;

    // Calculate the hypergeometric probability
    double p_value = (std::tgamma(n + 1) / (std::tgamma(a + 1) * std::tgamma(b + 1) * std::tgamma(c + 1) * std::tgamma(d + 1))) /
                     (std::tgamma(a + c + 1) / (std::tgamma(a + c) * std::tgamma(c + 1)) *
                      std::tgamma(b + d + 1) / (std::tgamma(b + d) * std::tgamma(d + 1)));

    // Note: This implementation only calculates the exact p-value for the given table.
    // For two-tailed test, one would need to calculate the probabilities for all tables
    // that are as extreme or more extreme than the observed one.

    return p_value;
}
