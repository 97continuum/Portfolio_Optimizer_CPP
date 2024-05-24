//
// Created by Talha Jamal on 20/05/2024.
//
#include "linearAlgebra.h"
#include <vector>

// Print Matrix
// Print Covariance Matrix
void printMatrix(const std::vector< std::vector<double> >& matrix, const std::string& matrixName)
{
    try
    {
        if (matrix.empty())
        {
            std::cerr << "Matrix is empty." << std::endl;
            return;
        }
        std::cout << "Matrix: " << matrixName << std::endl;
        for (const auto& row : matrix)
        {
            if (row.empty())
            {
                std::cerr << "Matrix " << matrixName << "contains an empty row." << std::endl;
                return;
            }
            for (const auto& value : row)
            {
                std::cout << value << " ";
            }
            std::cout << std::endl;
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Exception caught: " << e.what() << std::endl;
    }
}

// Print Vector
void printVector(const std::vector< double >& vector, const std::string& vectorName)
{
    try
    {
        if (vector.empty())
        {
            std::cerr << "Vector" << vectorName << "is empty" << std::endl;
            return;
        }
        std::cout << "Vector " << vectorName << std::endl;
        for (const auto& value : vector)
        {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Except caught: " << e.what() << std::endl;
    }
}

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

// Scalar Division
std::vector<double> scalarDivision(const std::vector<double>& vector, double scalar)
{
    std::vector<double> result(vector.size());
    for (size_t i = 0; i < vector.size(); ++i)
    {
        result[i] = vector[i] / scalar;
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

// Later on if Matrix Matrix addition/subtraction/division is needed, add that


// Unit Tests for each function

void testVectorAdd() {
    std::vector<double> a = {1.0, 2.0, 3.0};
    std::vector<double> b = {4.0, 5.0, 6.0};
    std::vector<double> result = vectorAddition(a, b);

    std::cout << "Vector Addition Result: ";
    for (double val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void testVectorSub() {
    std::vector<double> a = {4.0, 5.0, 6.0};
    std::vector<double> b = {1.0, 2.0, 3.0};
    std::vector<double> result = vectorSubtraction(a, b);

    std::cout << "Vector Subtraction Result: ";
    for (double val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void testDotProduct() {
    std::vector<double> a = {1.0, 2.0, 3.0};
    std::vector<double> b = {4.0, 5.0, 6.0};
    double result = vectorDotProduct(a, b);

    std::cout << "Dot Product Result: " << result << std::endl;
}

void testScalarMult() {
    std::vector<double> a = {1.0, 2.0, 3.0};
    double scalar = 2.0;
    std::vector<double> result = scalarMultiplication(a, scalar);

    std::cout << "Scalar Multiplication Result: ";
    for (double val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void testMatVecMult() {
    std::vector<std::vector<double>> mat = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0}
    };
    std::vector<double> vec = {1.0, 2.0, 3.0};
    std::vector<double> result = matrixVectorMultiplication(mat, vec);

    std::cout << "Matrix-Vector Multiplication Result: ";
    for (double val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

void testMatMatMult() {
    std::vector<std::vector<double>> A = {
            {1.0, 2.0},
            {3.0, 4.0},
            {5.0, 6.0}
    };
    std::vector<std::vector<double>> B = {
            {7.0, 8.0, 9.0},
            {10.0, 11.0, 12.0}
    };
    std::vector<std::vector<double>> result = matrixMultiplication(A, B);

    std::cout << "Matrix-Matrix Multiplication Result: " << std::endl;
    for (const auto& row : result) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

void testAllFunctions(){
    std::cout << "Testing Vector Addition..." << std::endl;
    testVectorAdd();
    std::cout << std::endl;

    std::cout << "Testing Vector Subtraction..." << std::endl;
    testVectorSub();
    std::cout << std::endl;

    std::cout << "Testing Dot Product..." << std::endl;
    testDotProduct();
    std::cout << std::endl;

    std::cout << "Testing Scalar Multiplication..." << std::endl;
    testScalarMult();
    std::cout << std::endl;

    std::cout << "Testing Matrix-Vector Multiplication..." << std::endl;
    testMatVecMult();
    std::cout << std::endl;

    std::cout << "Testing Matrix-Matrix Multiplication..." << std::endl;
    testMatMatMult();
    std::cout << std::endl;
}