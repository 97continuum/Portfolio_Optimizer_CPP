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
    int numberAssets = 5; // Initialize Number of Assets
    int numberReturns = 10; // Max Length of Returns Data

    string desiredDirectory = "/Users/talhajamal/Desktop/Code/Portfolio_Optimizer_CPP"; // Home Directory
    changeWorkingDirectory(desiredDirectory); // Change to Correct Working Directory
    string fileName = "data/asset_returns_small.csv"; // Return Data FileName
    checkFileInCurrentDirectory(fileName); // Check if File Exists and File Path is correct

    Matrix returnMatrix;
    returnMatrix.resize(numberReturns);
    for(int i=0;i<numberReturns;i++)
        returnMatrix[i].resize(numberAssets);
    readData(returnMatrix,fileName); // Read return data from the file and store in 2D returnMatrix

    Vector average = calculateAverage(returnMatrix);
    printVector(average, "average");
    Matrix covMatrix = calculateCovMatrix(returnMatrix);
    printMatrix(covMatrix, "Cov Matrix");

    // Set list of increasing target Portfolio Returns
    const int steps = 20;
    double temp[steps];
    double start = 0.005;
    double end = 0.1;
    double increment = (end - start) / (steps - 1);

    for (int i = 0; i < steps; ++i) {temp[i] = start + i * increment;}
    Vector tReturns(temp, temp + steps);

    ofstream resultsFile; // Create CSV File to send results to
    resultsFile.open("data/results.csv");
    resultsFile << "Target Returns,";
    for (int i = 0; i < numberAssets - 1; i++)
    {
        resultsFile << "Asset " << i+1 << ",";
    }
    resultsFile << numberAssets << endl;

    // Backtesting Parameters
    double epsilon = 10e-6;
    int isWindow = 100;
    int oosWindow = 12;
    int slidingWindow = 12;
    int numOfSlidingWindows = (numberReturns - isWindow - oosWindow) / (slidingWindow);


    /* Outer loop for each backtesting range
        // In Sample Start and End Date for Portfolio Building
        // Generate 20 Portfolios
        // Calculate return of each portfolio
        // Find average return by averaging the return of each portfolio
        // calculate actual average return: rT*w
        // Find covariance matrix for OOS data
    */

    int startIS, endIS, startOOS, endOOS; // indexes to control sliding window

    /*
    // Outer loop for each backtesting range
    for (int i = 0; i <= numOfSlidingWindows; i++)
    {
        startIS = (slidingWindow * i);
        endIS = startIS + (isWindow - 1);
        startOOS = startIS + isWindow;
        endOOS = startOOS + (oosWindow - 1);
        //cout << "Start IS: " << startIS << " End IS: " << endIS << " Start OOS: " << startOOS << " End OOS: " << endOOS;
        //cout << endl;
        // Run Portfolio Optimizer over every Return Target
        for (int j = 0; j < 1; j++) // change to 20 later
        {
            cout << "For Portfolio with Target Return " << tReturns[j]*100 << "%" << endl;
            cout << "Start IS: " << startIS << " End IS: " << endIS << endl;
            Matrix inSampleReturns = sliceMatrixByRows(returns, startIS, endIS);

            //Portfolio portfolio(inSampleReturns, tReturns[j], numberAssets, );
            // Design
            /* Initialize a Portfolio Object with spe
             *
             *
        }

    } */

/*
    Portfolio portfolio(returns, targetReturns, numberAssets, numberReturns);
    Vector portfolioMeanReturns = portfolio.calculateMeanReturn(); // Calculate mean returns
    //printVector(portfolioMeanReturns, "Mean Returns");
    Matrix portfolioCovMatrix = portfolio.calculateCovarianceMatrix(); // Calculate Cov Matrix
    //printMatrix(portfolioCovMatrix, "Cov Matrix"); // Print Covariance Matrix

    Vector weights = portfolio.solveOptimization();
    printVector(weights, "of Portfolio Weights - last two entries are Lagrange multiplier constants");
*/
    return 0;
}
