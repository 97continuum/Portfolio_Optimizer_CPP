//
// Created by Talha Jamal on 20/05/2024.
//

#ifndef COURSEWORK_PORTFOLIO_H
#define COURSEWORK_PORTFOLIO_H

#include <iostream>
#include "linearAlgebra.h"

using namespace std;

Matrix sliceMatrixByRows(const Matrix& matrix, int row_start, int row_end);
Vector calculateAverage(const Matrix& m);
Matrix calculateCovMatrix(const Matrix& m, Vector meanReturns);

class BacktestingEngine {
public:
    BacktestingEngine(int isWindow_,int oosWindow_,int slidingWindow_,int numOfSlidingWindows_,
                      const Matrix& AssetReturns_, double targetReturn_);

    // OOS Functions
    void calculateOOSMean();
    void calculateOOSCovMatrix();

    // IS Functions
    void calculateIsMean();
    void calculateIsCovMat();
    void calculateIsQ();
    void optimizer(double epsilon);

    void setb(double targetReturn_)
    {
        targetReturn = targetReturn_;
        Vector b(numAssets, 0.0);
        b.push_back(-targetReturn);
        b.push_back(-1);
        isb = b;
    }
    ~BacktestingEngine();
private:
    int isWindow, oosWindow, slidingWindow, numOfSlidingWindows;
    Matrix AssetReturns;
    size_t numAssets, numReturns;

    // OOS Variables
    Vector oosMean;
    vector<Matrix> oosCovMatrix;

    // IS Variables - including variables needed for Portfolio Optimization
    double targetReturn;
    Vector isb, isLambda, isMu;
    Matrix isMean, isWeights, X;
    vector<Matrix> isCovMatrix, Q;
};


class Portfolios : public BacktestingEngine {
public:
    Portfolios(int isWindow_, int oosWindow_, int slidingWindow_, int numOfSlidingWindows_,
               const Matrix& AssetReturns_, double targetReturn_);

private:
    Vector actualAverageReturn;
    Matrix actualCovMatrix;
    double annualizedActualReturn;
    double cumulativeActualReturn;
    double standardDeviation;
    double portfolioSharpeRatio;
};


/*
class Portfolio {
public:
    Portfolio(const vector< Vector >& returns, double targetReturn, int numAssets, int numPeriods);

    Vector calculateMeanReturn(); // Function to calculate
    Matrix calculateCovarianceMatrix(); // function to calculate Covariance Matrix
    Vector solveOptimization(); // Function to create Qx = b systems of linear equations
    static Vector conjugateGradient(const Matrix &Q, const Vector &b, const Vector &x0) ; // Function to solve Qx = b equation

private:
    Matrix returns; // Matrix to Store Returns
    double targetReturn; // Variable to hold Target Return for Portfolio
    int numAssets; // Total Number of Assets
    int numPeriods; // Total Number of Periods
    Vector meanReturns; // The vector of Mean Returns for N Assets
    Matrix covarianceMatrix; // Variance Covariance Matrix for N Assets
};*/

#endif //COURSEWORK_PORTFOLIO_H
