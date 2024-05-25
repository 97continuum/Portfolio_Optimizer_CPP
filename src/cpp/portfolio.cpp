//
// Created by Talha Jamal on 20/05/2024.
//

#include "portfolio.h"
#include <numeric>
#include "linearAlgebra.h"
#include <stdexcept>

using namespace std;

// Constructor for Portfolio Class
Portfolio::Portfolio(const Matrix &returns, double targetReturn, int numAssets, int numPeriods)
    : returns(returns), targetReturn(targetReturn), numAssets(numAssets), numPeriods(numPeriods){}

// Calculate Mean Returns
Vector Portfolio::calculateMeanReturn()
{
    //size_t numAssets = returns.size();
    //size_t numPeriods = returns[0].size();
    meanReturns.resize(numAssets);

    for (int i = 0; i < numAssets; ++i)
    {
        meanReturns[i] = accumulate(returns[i].begin(), returns[i].end(), 0.0) / numPeriods;
    }
    return meanReturns;
}

// Calculate Covariance Matrix
MatrixPortfolio::calculateCovarianceMatrix() {
    //size_t numAssets = returns.size();
    //size_t numPeriods = returns[0].size();
    covarianceMatrix.resize(numAssets, Vector(numAssets, 0.0));

    // Calculate Mean Returns if they have not been calculated
    if (meanReturns.empty()) {
        calculateMeanReturn();
    }

    for (size_t i = 0; i < numAssets; ++i) {
        for (size_t j = 0; j < numAssets; ++j) {
            double cov = 0.0;
            for (size_t k = 0; k < numPeriods; ++k) {
                cov += (returns[i][k] - meanReturns[i]) * (returns[j][k] - meanReturns[j]);
            }
            covarianceMatrix[i][j] = cov / (numPeriods - 1);
        }
    }
    return covarianceMatrix;
}

// solve Optimization Problem by creating the System of Linear Equations needed for Conjugate Gradient Method
Vector Portfolio::solveOptimization()
{
    //size_t n = meanReturns.size(); // number of Assets
    MatrixQ(numAssets + 2, Vector(numAssets + 2, 0.0)); // Dimensions of Matrix Q
    Vector b(numAssets + 2, 0.0); // Dimensions of Vector b
    Vector x0(numAssets, 1.0/numAssets); // Dimensions of Vector X
    x0.push_back(0.005);
    x0.push_back(0.005);

    // Fill Q Matrix
    for (size_t i=0; i < numAssets; ++i)
    {
        for (size_t j = 0; j < numAssets; ++j)
        {
            Q[i][j] = covarianceMatrix[i][j];
        }
        Q[i][numAssets] = -meanReturns[i];
        Q[i][numAssets+1] = -1.0;
        Q[numAssets][i] = -meanReturns[i];
        Q[numAssets+1][i] = -1.0;
    }
    // Fill b vector
    b[numAssets] = -targetReturn;
    b[numAssets+1] = -1.0;
    //printMatrix(Q, "Q");
    //printVector(b, "b");

    // Solve for Qx = b by Conjugate Method
    return conjugateGradient(Q, b, x0);
}

Vector Portfolio::conjugateGradient(const Matrix &Q,
                                                 const Vector &b,
                                                 const Vector &x0)
{
    Vector  x_, s_pre, s_aft, p_pre, p_aft, Qp, alphaQp, weights;
    double alpha, beta, pTQp;
    size_t n = b.size();

    for (size_t k = 0; k < n; ++k)
    {
        x_ = x0;
        s_pre = vectorSubtraction(b, matrixVectorMultiplication(Q, x0)); // s: b - Qx
        p_pre = s_pre; // set Initial Direction
        //sTs = vectorDotProduct(s_pre, s_pre);
        //printVector(x, "x");
        //printVector(s, "s");
        //printVector(p, "p");
        //cout << "sTs : " << sTs << endl;

        while (vectorDotProduct(s_pre, s_pre) > 1.0E-6) {
            Qp = matrixVectorMultiplication(Q, p_pre); //printVector(Qp, "Q*p");
            pTQp = vectorDotProduct(p_pre, Qp); //cout << "Dot Product " << pTQp << endl;
            alpha = vectorDotProduct(s_pre, s_pre) / pTQp; // step size //cout << "alpha " << alpha << endl;
            x_ = vectorAddition(x_, scalarMultiplication(p_pre, alpha)); //printVector(x, "x");
            alphaQp = scalarMultiplication(Qp, alpha);
            s_aft = vectorSubtraction(s_pre, alphaQp); //cout << "s_aft " << s_aft << endl;
            beta = vectorDotProduct(s_aft, s_aft) / vectorDotProduct(s_pre, s_pre);
            p_aft = vectorAddition(s_aft, scalarMultiplication(p_pre, beta));
            s_pre = s_aft;
            p_pre = p_aft;
        }
        for (int j = 0; j < x_.size()-2; j++)
        {
            weights.push_back(x_[j]); // Add weights except Langrange Multipliers from x
        }
    }
    //cout << "Finished Optimization" << endl;
    return weights;
}
// Helper Functions

// Slice Matrix
template <typename T>
vector< vector<T> > sliceMatrixByRows(const vector<vector<int>>& matrix,
                                                int row_start, int row_end) {
    // Check for valid range
    if (row_start > row_end || row_end > matrix.size())
    {
        throw out_of_range("Invalid slice range");
    }
    vector<vector<int>> submatrix(matrix.begin() + row_start, matrix.begin() + row_end);
    return submatrix;
}


// Backtesting Class