#ifndef QUANTITATIVE_ANALYSIS_HPP
#define QUANTITATIVE_ANALYSIS_HPP

#include <Eigen/Dense>
#include <vector>
#include <cmath>   // for sqrt and erf
#include <iostream>
#include <string>
#include <map>
#include <numeric>

// Function to perform linear regression and return beta coefficients, standard errors, and p-values
void linearRegression(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Eigen::VectorXd& beta, Eigen::VectorXd& se, Eigen::VectorXd& p_values);

// Function to calculate p-values based on the t-statistic
double calculatePValue(double t_stat, int degrees_of_freedom);

std::vector<std::vector<int>> create_binary_table(const std::unordered_map<std::string, bool>& groups, const std::vector<std::string>& list_path_snarl);

#endif // LINEAR_REGRESSION_HPP
