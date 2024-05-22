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
#include "parameter_estimation.h"
#include "linearAlgebra.h"
#include "unitTests.h"
#include "portfolio.h"

int  main (int  argc, char  *argv[])
{
    int numberAssets = 83; // Initialize Number of Assets
    int numberReturns = 700; // Max Length of Returns Data
    double **returnMatrix = new double*[numberAssets]; // a matrix to store the return data

    //allocate memory for return data
    for(int i=0;i<numberAssets;i++)
        returnMatrix[i]=new double[numberReturns];

    cout << "Reading Data" << std::endl;

    //read the data from the file and store it into the return matrix
    std::string fileName = "asset_returns.csv";
    checkFileInCurrentDirectory(fileName); // Check if File Exists and File Path is correct
    readData(returnMatrix,fileName); // Read return data from the file and store in 2D returnMatrix

    // Convert to vector of vectors
    std::vector<std::vector<double>> returns = convertToVectorMatrix(returnMatrix, numberAssets, numberReturns);

    // Test Linear Algebra Functions
    //testAllFunctions();

    std::cout << "Checking if both meanReturn and covariance matrix function are calculated correctly" << std::endl;

    std::cout << "Parameter Estimation Script:";
    std::cout << " =========================== ";
    std::cout << "Mean Return Calculation for First 10 Assets";
    int inSampleSize = 700;
    std::vector<double> meanReturns = calculateMean(returnMatrix, 10, inSampleSize);
    std::cout << "Mean Returns: " << std::endl;
    for (double mean: meanReturns)
        {
            std::cout << mean << " ,";
        }
    std::cout << std::endl;

    std::cout << "Portfolio Class Method:";
    std::cout << " =========================== ";
    std::cout << "Mean Return Calculation for First 10 Assets";

    double targetReturns = 0.10;
    Portfolio portfolio(returns, targetReturns, 10, inSampleSize);

    // Calculate mean returns
    std::vector<double> portfolioMeanReturns = portfolio.calculateMeanReturn();
    std::cout << "Portfolio Mean Returns:\n";
    for (double r : portfolioMeanReturns) {
        std::cout << r << " ";
    }
    std::cout << std::endl;


    std::cout << "Parameter Estimation Script:";
    std::cout << " =========================== ";
    std::cout << "Covariance Matrix Calculation for First 10 Assets";
    std::vector< std::vector<double> > covarianceMatrix = calculateCovarianceMatrix(returnMatrix, 10, inSampleSize);
    std::cout << "CovarianceMatrix" << std::endl;
    for (const auto& row: covarianceMatrix)
    {
        for(const auto& value: row)
        {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Portfolio Class Method:";
    std::cout << " =========================== ";
    std::cout << "Covariance Matrix Calculation for First 10 Assets";
    std::vector< std::vector<double> > portfolioCovMatrix = portfolio.calculateCovarianceMatrix();
    std::cout << "Covariance Matrix:" << std::endl;
    for (const auto& row: portfolioCovMatrix)
    {
        for (const auto& value: row)
        {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    // Delete Memory
    for(int i=0;i<numberAssets;i++)
        delete[] returnMatrix[i];
    delete[] returnMatrix;
    cout << "Deleted Memory" << std::endl;
    return 0;
}
