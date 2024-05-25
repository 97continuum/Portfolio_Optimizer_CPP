//
// Created by Talha Jamal on 20/05/2024.
//

#ifndef COURSEWORK_LINEARALGEBRA_H
#define COURSEWORK_LINEARALGEBRA_H

#include <iostream>
#include "vector"

using namespace std;

typedef vector<double> Vector;
typedef vector<Vector> Matrix;

void printMatrix(const Matrix& matrix, const string& matrixName);
void printVector(const Vector& vector, const string& vectorName);
Vector vectorAddition(const Vector& a, const Vector& b);
Vector vectorSubtraction(const Vector& a, const Vector& b);
double vectorDotProduct(const Vector& a, const Vector& b);
Vector scalarMultiplication(const Vector& vector, double scalar);
Vector scalarDivision(const Vector& vector, double scalar);
Vector matrixVectorMultiplication(const Matrix& matrix, const Vector& vector);
double vectorNorm(const Vector& vector);
Matrix matrixMultiplication(const Matrix& A, const Matrix& B);

// unit tests
void testVectorAdd();
void testVectorSub();
void testDotProduct();
void testScalarMult();
void testMatVecMult();
void testMatMatMult();
void testAllFunctions();

#endif //COURSEWORK_LINEARALGEBRA_H
