//
// Created by Talha Jamal on 20/05/2024.
//

#ifndef COURSEWORK_LINEARALGEBRA_H
#define COURSEWORK_LINEARALGEBRA_H

#include <iostream>
#include "vector"

void printMatrix(const std::vector< std::vector<double> >& matrix, const std::string& matrixName);
void printVector(const std::vector< double >& vector, const std::string& vectorName);
std::vector<double> vectorAddition(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> vectorSubtraction(const std::vector<double>& a, const std::vector<double>& b);
double vectorDotProduct(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> scalarMultiplication(const std::vector<double>& vector, double scalar);
std::vector<double> scalarDivision(const std::vector<double>& vector, double scalar);
std::vector<double> matrixVectorMultiplication(const std::vector< std::vector<double> >& matrix,
                                               const std::vector<double>& vector);
double vectorNorm(const std::vector<double>& vector);
std::vector< std::vector<double> > matrixMultiplication(const std::vector< std::vector<double> >& A,
                                                      const std::vector< std::vector<double> >& B);

// unit tests
void testVectorAdd();
void testVectorSub();
void testDotProduct();
void testScalarMult();
void testMatVecMult();
void testMatMatMult();
void testAllFunctions();

#endif //COURSEWORK_LINEARALGEBRA_H
