#include "Utilities.h"

// Square function
double utl::sqr(const double &number)
{
    return number * number;
}

// Vector of ones
vecDbl utl::ones(const size_t &n)
{
    return vecDbl(n, 1.0);
}

Matrix utl::ones(const size_t &nrow, const size_t &ncol)
{
    Matrix mat;
    mat.reserve(nrow);
    for (size_t i = 0u; i < nrow; ++i)
        mat.push_back(utl::ones(ncol));
    return mat;
}

// Vector of zeros
vecDbl utl::zeros(const size_t &n)
{
    return vecDbl(n, 0.0);
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
    return vecUns(n, 0u);
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
    return vecDbl(n, x);
}

// Repeat an unsigned integer many times
vecUns utl::repUns(const size_t &x, const size_t &n)
{
    return vecUns(n, x);
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
    for (size_t i = 0u; i < nrow; ++i) assert(m[i].size() == ncol);

    // Check that the last row and the last column are full of zeros
    for (size_t i = 0u; i < nrow; ++i) assert(m[i][ncol - 1u] == 0.0);
    for (size_t j = 0u; j < ncol; ++j) assert(m[nrow - 1u][j] == 0.0);

    // Calculate row sums
    for (size_t i = 0u; i < nrow; ++i) {
        double sum = 0.0;
        for (size_t j = 0u; j < ncol - 1u; ++j) {
            sum += m[i][j];
            m[i][ncol - 1u] = sum;
        }
    }

    // Calculate column sums
    for (size_t j = 0u; j < ncol; ++j) {
        double sum = 0.0;
        for (size_t i = 0u; i < nrow - 1u; ++i) {
            sum += m[i][j];
            m[nrow - 1u][j] = sum;
        }
    }
}

// Same for a matrix of unsigned integers
void utl::marginalize(MatUns &m)
{
    const size_t nrow = m.size();
    const size_t ncol = m[0u].size();

    // Check that all rows have the same number of elements
    for (size_t i = 0u; i < nrow; ++i) assert(m[i].size() == ncol);

    // Check that the last row and the last column are full of zeros
    for (size_t i = 0u; i < nrow; ++i) assert(m[i][ncol - 1u] == 0u);
    for (size_t j = 0u; j < ncol; ++j) assert(m[nrow - 1u][j] == 0u);

    // Calculate row sums
    for (size_t i = 0u; i < nrow; ++i) {
        size_t sum = 0u;
        for (size_t j = 0u; j < ncol - 1u; ++j) {
            sum += m[i][j];
            m[i][ncol - 1u] = sum;
        }
    }

    // Calculate column sums
    for (size_t j = 0u; j < ncol; ++j) {
        size_t sum = 0u;
        for (size_t i = 0u; i < nrow - 1u; ++i) {
            sum += m[i][j];
            m[nrow - 1u][j] = sum;
        }
    }
}

// Pairwise division between two matrices (the second contains unsigned int)
Matrix utl::dividemat(const Matrix &m1, const MatUns &m2)
{

    const size_t nrow = m1.size();
    const size_t ncol = m1[0u].size();

    // Check dimensions
    assert(nrow == m2.size());
    assert(ncol == m1[0u].size());
    for (size_t i = 0u; i < nrow; ++i) {
        assert(m1[i].size() == ncol);
        assert(m2[i].size() == m2[0u].size());
    }

    // Make an empty matrix
    Matrix res = utl::zeros(nrow, ncol);

    // Fill it in with pairwise ratios
    for (size_t i = 0u; i < nrow; ++i) {
        for (size_t j = 0u; j < ncol; ++j) {
            if (m2[i][j] != 0u) {
                res[i][j] = m1[i][j] / m2[i][j];
            }
        }
    }

    return res;
}

// Round to the ith decimal
double utl::round(const double &x, const size_t &i)
{
    double p = 1.0;
    for (size_t j = 0u; j < i; ++j) p *= 10.0;
    return std::round(x * p) / p;
}

void utl::correct(double &x, const double &x0, const double &d)
{
    x = x < x0 + d && x > x0 - d ? x0 : x;
}

// Convert unsigned integer to double
double utl::size2dbl(const size_t &x)
{
    return static_cast<double>(x);
}

// Convert double to unsigned integer
size_t utl::dbl2size(const double &x)
{
    return static_cast<size_t>(x);
}

// Save to file
void stf::write(const double &x, std::shared_ptr<std::ofstream> &out)
{
    out->write((char *) &x, sizeof(x));
}

void stf::write(const vecDbl &v, std::shared_ptr<std::ofstream> &out)
{
    if (v.size() > 0u)
        for (auto x : v)
            stf::write(x, out);
}
