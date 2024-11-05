#pragma once

#ifndef BINARY_ANALYSIS_HPP
#define BINARY_ANALYSIS_HPP

#include <Eigen/Dense>
#include <vector>

// Function to perform the Chi-Squared Test
// Returns the Chi-Squared statistic and sets the p-value reference
double chiSquaredTest(const Eigen::MatrixXd& observed, double& p_value);

// Function to perform Fisher's Exact Test
// Returns the p-value for the Fisher's Exact Test
double fisherExactTest(const Eigen::MatrixXd& contingencyTable);

#endif
