//
// Created by Talha Jamal on 19/05/2024.
//
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include "csv.h"
#include "read_data.h"
#include "linearAlgebra.h"
#include "portfolio.h"

using namespace std;

int  main (int  argc, char  *argv[])
{
    int numberAssets = 83; // Initialize Number of Assets
    int numberReturns = 700; // Max Length of Returns Data

    string desiredDirectory = "/Users/talhajamal/Desktop/Code/Portfolio_Optimizer_CPP"; // Home Directory
    changeWorkingDirectory(desiredDirectory); // Change to Correct Working Directory
    string fileName = "data/asset_returns.csv"; // Return Data FileName
    checkFileInCurrentDirectory(fileName); // Check if File Exists and File Path is correct

    Matrix returnMatrix;
    returnMatrix.resize(numberReturns);
    for(int i=0;i<numberReturns;i++)
        returnMatrix[i].resize(numberAssets);
    readData(returnMatrix,fileName); // Read return data from the file and store in 2D returnMatrix

    // Vector average = calculateAverage(returnMatrix);
    // printVector(average, "average");
    // Matrix covMatrix = calculateCovMatrix(returnMatrix, average);
    // printMatrix(covMatrix, "Cov Matrix");

    // Set list of increasing target Portfolio Returns
    const int steps = 20;
    double temp[steps];
    double start = 0.005;
    double end = 0.1;
    double increment = (end - start) / (steps - 1);

    for (int i = 0; i < steps; ++i) {temp[i] = start + i * increment;}
    Vector tReturns(temp, temp + steps);

    // Backtesting Parameters
    double epsilon = 10e-6;
    int isWindow = 100;
    int oosWindow = 12;
    int slidingWindow = 12;
    int numOfSlidingWindows = (numberReturns - isWindow - oosWindow) / (slidingWindow);

    ofstream resultsFile; // Create CSV File to send results to
    resultsFile.open("data/results.csv");
    resultsFile << "Target Returns,";
    for (int i = 0; i < numOfSlidingWindows - 1; ++i)
    {
        resultsFile << "Backtest Period: " << i+1 << ",";
    }
    resultsFile << "Backtest Period: " << numOfSlidingWindows << endl;

    vector<Portfolios> TargetReturnsPortfolios;

    for (int i = 0; i < 20; i++)
    {
        cout << "Running a Backtest for Portfolio with Target Returns of : " << tReturns[i]*100 << endl;
        Portfolios portfolio(isWindow, oosWindow, slidingWindow, numOfSlidingWindows, returnMatrix, tReturns[i]);
        portfolio.calculateIsMean();
        portfolio.calculateIsCovMat();
        portfolio.calculateOOSMean();
        portfolio.calculateOOSCovMatrix();
        portfolio.calculateIsQ();
        portfolio.optimizer(epsilon);
        portfolio.runBacktest();
        TargetReturnsPortfolios.push_back(portfolio);

        cout << "With a Target Return of " << tReturns[i]*100 << " the actual return of the Portfolio is: " << endl;
        Vector actualAvgReturn = portfolio.getActualAverageReturn();
        printVector(actualAvgReturn, "Actual Average Returns");
        cout << "Portfolio's Sharpe Ratio is: " << portfolio.getPortfolioSharpeRatio() << endl;
        cout << "Portfolio's Actual Average Abnormal Return is: " << portfolio.getAvgAbnormalReturn() << endl;
        cout << "Portfolio's Cumulate Average Return is: " << portfolio.getCumulativeAvgAbnormalReturn() << endl;

        resultsFile << tReturns[i] << ",";
        for (int j = 0; j < actualAvgReturn.size()-1; j++){
            resultsFile << actualAvgReturn[j] << ",";
        }
        resultsFile << actualAvgReturn[actualAvgReturn.size()-1] << endl;
    }
    return 0;
}
