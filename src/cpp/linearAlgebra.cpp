//
// Created by Talha Jamal on 20/05/2024.
//
#include "linearAlgebra.h"
#include <vector>

using namespace std;

// Print Matrix
// Print Covariance Matrix
void printMatrix(const Matrix & matrix, const string& matrixName)
{
    try
    {
        if (matrix.empty())
        {
            cerr << "Matrix is empty." << endl;
            return;
        }
        cout << "Matrix: " << matrixName << endl;
        for (const auto& row : matrix)
        {
            if (row.empty())
            {
                cerr << "Matrix " << matrixName << "contains an empty row." << endl;
                return;
            }
            for (const auto& value : row)
            {
                cout << value << " ";
            }
            cout << endl;
        }
    }
    catch (const exception& e)
    {
        cerr << "Exception caught: " << e.what() << endl;
    }
}

// Print Vector
void printVector(const Vector& vector, const string& vectorName)
{
    try
    {
        if (vector.empty())
        {
            cerr << "Vector" << vectorName << "is empty" << endl;
            return;
        }
        cout << "Vector " << vectorName << endl;
        for (const auto& value : vector)
        {
            cout << value << " ";
        }
        cout << endl;
    }
    catch (const exception& e)
    {
        cerr << "Except caught: " << e.what() << endl;
    }
}

// Vector Addition
Vector vectorAddition(const Vector& a, const Vector& b)
{
    Vector result(a.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        result[i] = a[i] + b[i];
    }
    return result;
}

// Vector Subtraction
Vector vectorSubtraction(const Vector& a, const Vector& b)
{
    Vector result(a.size());
    for (size_t i = 0; i < a.size(); ++i)
    {
        result[i] = a[i] - b[i];
    }
    return result;
}

// Vector Dot Product
double vectorDotProduct(const Vector& a, const Vector& b)
{
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        result += a[i] * b[i];
    }
    return result;
}

// Vector Scalar Multiplication
Vector scalarMultiplication(const Vector& vector, double scalar)
{
    Vector result(vector.size());
    for (size_t i = 0; i < vector.size(); ++i)
    {
        result[i] = scalar * vector[i];
    }
    return result;
}

// Vector Scalar Division
Vector scalarDivision(const Vector& vector, double scalar)
{
    Vector result(vector.size());
    for (size_t i = 0; i < vector.size(); ++i)
    {
        result[i] = vector[i] / scalar;
    }
    return result;
}

// Vector Norm
double vectorNorm(const Vector& vector)
{
    return sqrt(vectorDotProduct(vector, vector));
}

// Matrix Vector Multiplication
Vector matrixVectorMultiplication(const Matrix& matrix, const Vector& vector)
{
    // Check if the matrix is empty
    if (matrix.empty() || vector.empty()) {
        throw invalid_argument("Matrix or vector is empty.");
    }
    // Check if the number of columns in the matrix matches the size of the vector
    if (matrix[0].size() != vector.size()) {
        throw invalid_argument("Matrix columns do not match vector size.");
    }
    Vector result(matrix.size(), 0.0);  // Initialize the result vector with zeros

    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < vector.size(); ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}


Matrix matrixMultiplication(const Matrix& A, const Matrix& B)
{
    size_t n = A.size(); // Number of rows in A
    size_t m = A[0].size(); // Number of Columns in A
    size_t p = B[0].size(); // Number of Columns in B

    // Check if the Number of columns of A match the number of Rows in B
    if (m != B.size())
    {
        throw invalid_argument("Dimensions of Matrix A Columns do not match Dimensions of Matrix B Rows");
    }

    Matrix result(n, Vector(p, 0.0)); // Initialize result matrix with zeros
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
    Vector a = {1.0, 2.0, 3.0};
    Vector b = {4.0, 5.0, 6.0};
    Vector result = vectorAddition(a, b);

    cout << "Vector Addition Result: ";
    for (double val : result) {
        cout << val << " ";
    }
    cout << endl;
}

void testVectorSub() {
    Vector a = {4.0, 5.0, 6.0};
    Vector b = {1.0, 2.0, 3.0};
    Vector result = vectorSubtraction(a, b);

    cout << "Vector Subtraction Result: ";
    for (double val : result) {
        cout << val << " ";
    }
    cout << endl;
}

void testDotProduct() {
    Vector a = {1.0, 2.0, 3.0};
    Vector b = {4.0, 5.0, 6.0};
    double result = vectorDotProduct(a, b);

    cout << "Dot Product Result: " << result << endl;
}

void testScalarMult() {
    Vector a = {1.0, 2.0, 3.0};
    double scalar = 2.0;
    Vector result = scalarMultiplication(a, scalar);

    cout << "Scalar Multiplication Result: ";
    for (double val : result) {
        cout << val << " ";
    }
    cout << endl;
}

void testMatVecMult() {
    Matrix mat = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0}
    };
    Vector vec = {1.0, 2.0, 3.0};
    Vector result = matrixVectorMultiplication(mat, vec);

    cout << "Matrix-Vector Multiplication Result: ";
    for (double val : result) {
        cout << val << " ";
    }
    cout << endl;
}

void testMatMatMult() {
    Matrix A = {
            {1.0, 2.0},
            {3.0, 4.0},
            {5.0, 6.0}
    };
    Matrix B = {
            {7.0, 8.0, 9.0},
            {10.0, 11.0, 12.0}
    };
    Matrix result = matrixMultiplication(A, B);

    cout << "Matrix-Matrix Multiplication Result: " << endl;
    for (const auto& row : result) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}

void testAllFunctions(){
    cout << "Testing Vector Addition..." << endl;
    testVectorAdd();
    cout << endl;

    cout << "Testing Vector Subtraction..." << endl;
    testVectorSub();
    cout << endl;

    cout << "Testing Dot Product..." << endl;
    testDotProduct();
    cout << endl;

    cout << "Testing Scalar Multiplication..." << endl;
    testScalarMult();
    cout << endl;

    cout << "Testing Matrix-Vector Multiplication..." << endl;
    testMatVecMult();
    cout << endl;

    cout << "Testing Matrix-Matrix Multiplication..." << endl;
    testMatMatMult();
    cout << endl;
}