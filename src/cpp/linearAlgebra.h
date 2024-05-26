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

Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator*(const Vector& vector, double scalar);
Vector operator/(const Vector& vector, double scalar);
Vector operator*(const Matrix& matrix, const Vector& vector);
double operator*(const Vector& a, const Vector& b);

Matrix operator+(const Matrix& a, const Matrix& b);
Matrix operator-(const Matrix& a, const Matrix& b);
Matrix operator*(const Matrix& a, double b);
Matrix operator*(const Matrix& A, const Matrix& B);
Matrix operator/(const Matrix& a, double b);

// unit tests
void testVectorAdd();
void testVectorSub();
void testDotProduct();
void testScalarMult();
void testMatVecMult();
void testMatMatMult();
void testAllFunctions();

void printMatrix(const Matrix& matrix, const string& matrixName);
void printVector(const Vector& vector, const string& vectorName);

#endif //COURSEWORK_LINEARALGEBRA_H
