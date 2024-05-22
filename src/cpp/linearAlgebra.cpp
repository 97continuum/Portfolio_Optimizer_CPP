//
// Created by Talha Jamal on 20/05/2024.
//
#include "linearAlgebra.h"
#include <vector>

// Vector Addition
std::vector<double> vectorAddition(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        result[i] = a[i] + b[i];
    }
    return result;
}

// Vector Subtraction
std::vector<double> vectorSubtraction(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        result[i] = a[i] - b[i];
    }
    return result;
}

// Dot Product
double vectorDotProduct(const std::vector<double>& a, const std::vector<double>& b)
{
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

// Scalar Multiplication
std::vector<double> scalarMultiplication(const std::vector<double>& vector, double scalar)
{
    std::vector<double> result(vector.size());
    for (size_t i = 0; i < vector.size(); ++i)
    {
        result[i] = scalar * vector[i];
    }
    return result;
}

// Matrix Vector Multiplication
std::vector<double> matrixVectorMultiplication(const std::vector< std::vector<double> >& matrix,
                                               const std::vector<double>& vector)
{
    // Check if the matrix is empty
    if (matrix.empty() || vector.empty()) {
        throw std::invalid_argument("Matrix or vector is empty.");
    }
    // Check if the number of columns in the matrix matches the size of the vector
    if (matrix[0].size() != vector.size()) {
        throw std::invalid_argument("Matrix columns do not match vector size.");
    }
    std::vector<double> result(matrix.size(), 0.0);  // Initialize the result vector with zeros

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < vector.size(); ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

double vectorNorm(const std::vector<double>& vector)
{
    return std::sqrt(vectorDotProduct(vector, vector));
}

std::vector< std::vector<double> > matrixMultiplication(const std::vector< std::vector<double> >& A,
                                                      const std::vector< std::vector<double> >& B)
{
    size_t n = A.size(); // Number of rows in A
    size_t m = A[0].size(); // Number of Columns in A
    size_t p = B[0].size(); // Number of Columns in B

    // Check if the Number of columns of A match the number of Rows in B
    if (m != B.size())
    {
        throw std::invalid_argument("Dimensions of Matrix A Columns do not match Dimensions of Matrix B Rows");
    }

    // Initialize result matrix with zeros
    std::vector< std::vector<double> > result(n, std::vector<double>(p, 0.0));

    // Perform matrix multiplication
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) {
            for (size_t k = 0; k < m; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}
