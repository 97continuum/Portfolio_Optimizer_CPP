//
// Created by Talha Jamal on 20/05/2024.
//

#ifndef COURSEWORK_PORTFOLIO_H
#define COURSEWORK_PORTFOLIO_H

#include <iostream>
#include "linearAlgebra.h"

using namespace std;

template <typename T>
vector< vector<T> > sliceMatrixByRows(const Matrix& matrix, int row_start, int row_end);

class Portfolio {
public:
    Portfolio(const vector< Vector >& returns, double targetReturn, int numAssets, int numPeriods);

    Vector calculateMeanReturn(); // Function to calculate
    MatrixcalculateCovarianceMatrix(); // function to calculate Covariance Matrix
    Vector solveOptimization(); // Function to create Qx = b systems of linear equations
    static Vector conjugateGradient(const Matrix &Q,
                                          const Vector &b,
                                          const Vector &x0) ; // Function to solve Qx = b equation

private:
    Matrixreturns; // Matrix to Store Returns
    double targetReturn; // Variable to hold Target Return for Portfolio
    int numAssets; // Total Number of Assets
    int numPeriods; // Total Number of Periods
    Vector meanReturns; // The vector of Mean Returns for N Assets
    MatrixcovarianceMatrix; // Variance Covariance Matrix for N Assets
};

#endif //COURSEWORK_PORTFOLIO_H
