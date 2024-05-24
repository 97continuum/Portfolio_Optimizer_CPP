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

int  main (int  argc, char  *argv[])
{
    std::string desiredDirectory = "/Users/talhajamal/Desktop/Code/Portfolio_Optimizer_CPP"; // Home Directory
    changeWorkingDirectory(desiredDirectory); // Change to Correct Working Directory

    // Initialize Variables
    int numberAssets = 83; // Initialize Number of Assets
    int numberReturns = 700; // Max Length of Returns Data
    auto **returnMatrix = new double*[numberAssets]; // a matrix to store the return data
    //allocate memory for return data
    for(int i=0;i<numberAssets;i++)
        returnMatrix[i]=new double[numberReturns];

    //read the data from the file and store it into the return matrix
    std::string fileName = "data/asset_returns.csv";
    checkFileInCurrentDirectory(fileName); // Check if File Exists and File Path is correct
    readData(returnMatrix,fileName); // Read return data from the file and store in 2D returnMatrix
    // Convert to vector of vectors
    std::vector<std::vector<double>> returns = convertToVectorMatrix(returnMatrix, numberAssets, numberReturns);

    // Set list of increasing target Portfolio Returns
    const int steps = 20;
    double temp[steps];
    double start = 0.005;
    double end = 0.1;
    double increment = (end - start) / (steps - 1);
    for (int i = 0; i < steps; ++i) {
        temp[i] = start + i * increment;
    }
    std::vector<double> tReturns(temp, temp + steps);
    for (double value : tReturns) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    // Backtesting Parameters
    int isWindow = 100;
    int oosWindow = 12;
    int slidingWindow = 12;
    int numOfSlidingWindows = (numberReturns - isWindow - oosWindow) / (slidingWindow);
    //std::cout << numOfSlidingWindows << std::endl;

    // Run backtest and send returns to CSV file
    ofstream resultsFile;
    resultsFile.open("data/results.csv");
    resultsFile << "Target Returns,";
    for (int i = 0; i < numberAssets - 1; i++){
        resultsFile << "Asset " << i+1 << ",";
    }
    resultsFile << numberAssets << endl;

    double epsilon = 10e-6;

    // Run Portfolio Optimizer over every Return Target
    for (int i = 0; i < 20; i++)
    {
        std::cout << "For Portfolio with Target Return " << tReturns[i]*100 << "%" << std::endl;

    }

/*
    Portfolio portfolio(returns, targetReturns, numberAssets, numberReturns);
    std::vector<double> portfolioMeanReturns = portfolio.calculateMeanReturn(); // Calculate mean returns
    //printVector(portfolioMeanReturns, "Mean Returns");
    std::vector< std::vector<double> > portfolioCovMatrix = portfolio.calculateCovarianceMatrix(); // Calculate Cov Matrix
    //printMatrix(portfolioCovMatrix, "Cov Matrix"); // Print Covariance Matrix

    std::vector<double> weights = portfolio.solveOptimization();
    printVector(weights, "of Portfolio Weights - last two entries are Langrage multiplier constants");
*/

    // Delete Memory from Double Pointer
    deleteDoublePointer(returnMatrix, numberAssets);
    return 0;
}
