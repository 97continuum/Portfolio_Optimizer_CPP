//
// Created by Talha Jamal on 20/05/2024.
//

#ifndef COURSEWORK_PARAMETER_ESTIMATION_H
#define COURSEWORK_PARAMETER_ESTIMATION_H

#include <vector>

std::vector<double> calculateMean(double** data, int numberAssets, int n);
std::vector< std::vector<double> > calculateCovarianceMatrix(double** data, int numberAssets, int n);

#endif //COURSEWORK_PARAMETER_ESTIMATION_H
