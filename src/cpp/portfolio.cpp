//
// Created by Talha Jamal on 20/05/2024.
//

#include <numeric>
#include "portfolio.h"
#include "linearAlgebra.h"
#include <stdexcept>
#include <cmath>

using namespace std;

// Slice Matrix
Matrix sliceMatrixByRows(const Matrix& matrix, int row_start, int row_end)
{
    cout << matrix[0].size() << endl;
    // Check for valid range
    if (row_start > row_end || row_end > matrix[0].size())
    {
        throw out_of_range("Invalid slice range");
    }
    Matrix submatrix(matrix.begin() + row_start, matrix.begin() + row_end);
    return submatrix;
}

Vector calculateAverage(const Matrix& m)
{
    if (m.empty()){return {};}
    size_t numRows = m.size();
    size_t numCols = m[0].size();
    Vector result(numCols, 0.0);
    // Sum each column's elements
    for (int i = 0; i < numRows; ++i)
    {
        for (int j = 0; j < numCols; ++j)
        {
            result[j] += m[i][j];
        }
    }
    // Divide by number of Rows to get the mean
    for (int i = 0; i < numCols; ++i)
    {
        result[i] /= numRows;
    }
    return result;
}

Matrix calculateCovMatrix(const Matrix& m, Vector meanReturns)
{
    size_t numCols = m[0].size();
    size_t numRows = m.size();
    Matrix covMatrix(numCols, Vector(numCols, 0.0));
    //Vector meanReturns = calculateAverage(m);

    for (size_t i = 0; i < numCols; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            double cov = 0.0;
            for (size_t k = 0; k < numRows; ++k) {
                cov += (m[k][i] - meanReturns[i]) * (m[k][j] - meanReturns[j]);
            }
            covMatrix[i][j] = cov / (numRows - 1);
        }
    }
    return covMatrix;
}

// Backtesting Class Constructor
BacktestingEngine::BacktestingEngine(int isWindow_,int oosWindow_,int slidingWindow_,int numOfSlidingWindows_,
                                     const Matrix& AssetReturns_, double targetReturn_):
                                     isWindow(isWindow_),oosWindow(oosWindow_),slidingWindow(slidingWindow_),
                                     numOfSlidingWindows(numOfSlidingWindows_),AssetReturns(AssetReturns_),
                                     targetReturn(targetReturn_)
{
    numReturns = AssetReturns.size();
    numAssets = AssetReturns[0].size();
    targetReturn = targetReturn_;
    Vector temp(numAssets, 0.0);
    temp.push_back(-targetReturn_);
    temp.push_back(-1);
    isb = temp;
}

// Destructor for Backtesting Engine Class
BacktestingEngine::~BacktestingEngine() = default;

// Constructor for Portfolios Class
Portfolios::Portfolios(int isWindow_, int oosWindow_, int slidingWindow_, int numOfSlidingWindows_,
                       const Matrix& AssetReturns_, double targetReturn_)
        : BacktestingEngine(isWindow_, oosWindow_, slidingWindow_, numOfSlidingWindows_, AssetReturns_, targetReturn_) {}




/*
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
Matrix Portfolio::calculateCovarianceMatrix() {
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
 */

/*
// solve Optimization Problem by creating the System of Linear Equations needed for Conjugate Gradient Method
Vector Portfolio::solveOptimization()
{
    //size_t n = meanReturns.size(); // number of Assets
    Matrix Q(numAssets + 2, Vector(numAssets + 2, 0.0)); // Dimensions of Matrix Q
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
    printMatrix(Q, "Q");
    printVector(x0, "x");
    printVector(b, "b");

    // Solve for Qx = b by Conjugate Method
    return conjugateGradient(Q, b, x0);
}
*/
/*
Vector Portfolio::conjugateGradient(const Matrix &Q, const Vector &b, const Vector &x0)
{
    Vector  x_, s_pre, s_aft, p_pre, p_aft, Qp, alphaQp, weights;
    double alpha, beta, pTQp;
    size_t n = b.size();

    for (size_t k = 0; k < n; ++k)
    {
        x_ = x0;
        s_pre = b - (Q * x0);// s: b - Qx
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
 */
// Helper Functions
