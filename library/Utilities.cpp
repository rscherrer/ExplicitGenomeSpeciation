#include "Utilities.h"

// Square function
double utl::sqr(const double &number)
{
    return number * number;
}

// Vector of ones
vecDbl utl::ones(const size_t &n)
{
    vecDbl ones;
    ones.reserve(n);
    for (size_t i = 0u; i < n; ++i)
        ones.push_back(1.0);
    return ones;
}

// Vector of zeros
vecDbl utl::zeros(const size_t &n)
{
    vecDbl zeros;
    zeros.reserve(n);
    for (size_t i = 0u; i < n; ++i)
        zeros.push_back(0.0);
    return zeros;
}

// Matrix of zeros
Matrix utl::zeros(const size_t &nrow, const size_t &ncol)
{
    Matrix mat;
    mat.reserve(nrow);
    for (size_t i = 0u; i < nrow; ++i)
        mat.push_back(utl::zeros(ncol));
    return mat;
}

// 3D matrix of zeros
Matx3d utl::zeros(const size_t &n0, const size_t &n1, const size_t &n2)
{
    Matx3d mat;
    mat.reserve(n0);
    for (size_t i = 0u; i < n0; ++i)
        mat.push_back(utl::zeros(n1, n2));
    return mat;
}

// Unsigned zeros
vecUns utl::uzeros(const size_t &n)
{
    vecUns zeros;
    zeros.reserve(n);
    for (size_t i = 0u; i < n; ++i)
        zeros.push_back(0u);
    return zeros;
}

// Matrix of unsigned zeros
MatUns utl::uzeros(const size_t &nrow, const size_t &ncol)
{
    MatUns mat;
    mat.reserve(nrow);
    for (size_t i = 0u; i < nrow; ++i)
        mat.push_back(utl::uzeros(ncol));
    return mat;
}

// Repeat a number many times
vecDbl utl::rep(const double &x, const size_t &n)
{
    vecDbl reps;
    reps.reserve(n);
    for (size_t i = 0u; i < n; ++i)
        reps.push_back(x);
    return reps;
}

// Repeat an unsigned integer many times
vecUns utl::repUns(const size_t &x, const size_t &n)
{
    vecUns reps;
    reps.reserve(n);
    for (size_t i = 0u; i < n; ++i)
        reps.push_back(x);
    return reps;
}


// Sum of doubles
double utl::sum(vecDbl &v)
{
    return std::accumulate(v.begin(), v.end(), 0);
}

// Sum of a matrix of doubles
double utl::sum(Matrix &m)
{
    double res = 0.0;
    for (size_t i = 0u; i < m.size(); ++i)
        res += utl::sum(m[i]);
    return res;
}

// Find the position of the minimum element in a vector of floats
size_t utl::argmin(vecDbl &v)
{
    return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
}

// Sum of unsigned integers
size_t utl::sumu(const vecUns &v) {
    size_t sum = 0u;
    for (size_t i = 0u; i < v.size(); ++i) sum += v[i];
    return sum;
}

// Calculate sums in the margins of a matrix
void utl::marginalize(Matrix &m)
{
    const size_t nrow = m.size();
    const size_t ncol = m[0u].size();

    // Check that all rows have the same number of elements
    for (size_t i = 0u; i < nrow; ++i) {
        assert(m[i].size() == ncol);
    }

    // Calculate row sums
    for (size_t i = 0u; i < nrow; ++i) {
        for (size_t j = 0u; j < ncol - 1u; ++j) {
            m[i][ncol - 1u] += m[i][j];
        }
    }

    // Calculate column sums
    for (size_t j = 0u; j < ncol; ++j) {
        for (size_t i = 0u; i < nrow - 1u; ++i) {
            m[nrow - 1u][j] += m[i][j];
        }
    }
}

// Same for a matrix of unsigned integers
void utl::marginalize(MatUns &m)
{
    const size_t nrow = m.size();
    const size_t ncol = m[0u].size();

    // Check that all rows have the same number of elements
    for (size_t i = 0u; i < nrow; ++i) {
        assert(m[i].size() == ncol);
    }

    // Calculate row sums
    for (size_t i = 0u; i < nrow; ++i) {
        for (size_t j = 0u; j < ncol - 1u; ++j) {
            m[i][ncol - 1u] += m[i][j];
        }
    }

    // Calculate column sums
    for (size_t j = 0u; j < ncol; ++j) {
        for (size_t i = 0u; i < nrow - 1u; ++i) {
            m[nrow - 1u][j] += m[i][j];
        }
    }
}

// Pairwise division between two matrices (the second contains unsigned int)
Matrix utl::dividemat(const Matrix &m1, const MatUns &m2)
{

    const size_t nrow1 = m1.size();
    const size_t nrow2 = m2.size();
    const size_t ncol1 = m1[0u].size();
    const size_t ncol2 = m2[0u].size();

    // Check dimensions
    assert(nrow1 == nrow2);
    assert(ncol1 == ncol2);
    for (size_t i = 0u; i < nrow1; ++i) {
        assert(m1[i].size() == ncol1);
        assert(m2[i].size() == ncol2);
    }

    // Make an empty matrix
    Matrix res = utl::zeros(nrow1, ncol1);

    // Fill it in with pairwise ratios
    for (size_t i = 0u; i < nrow1; ++i) {
        for (size_t j = 0u; j < ncol1; ++j) {
            if (m2[i][j] != 0u) {
                res[i][j] = m1[i][j] / m2[i][j];
            }
        }
    }

    return res;
}

// Convert unsigned integer to double
double utl::size2dbl(const size_t &x)
{
    return static_cast<double>(x);
}
