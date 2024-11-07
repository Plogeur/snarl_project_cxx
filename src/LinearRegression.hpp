#ifndef LINEAR_REGRESSION_HPP
#define LINEAR_REGRESSION_HPP

#include <Eigen/Dense>
#include <vector>
#include <cmath>   // for sqrt and erf
#include <iostream>

// Function to perform linear regression and return beta coefficients, standard errors, and p-values
void linearRegression(const Eigen::MatrixXd& X, const Eigen::VectorXd& y, Eigen::VectorXd& beta, Eigen::VectorXd& se, Eigen::VectorXd& p_values);

// Function to calculate p-values based on the t-statistic
double calculatePValue(double t_stat, int degrees_of_freedom);

#endif // LINEAR_REGRESSION_HPP
