#include "Utilities.h"

double utl::Xst(const double &v0, const double &v1, const double &v,
 const size_t &n0, const size_t &n1, const size_t &n,
  const double &tiny)
{
    if (v < tiny) return 0.0;

    const double xst = 1.0 - (n0 * v0 + n1 * v1) / (n * v);

    if (xst < tiny) return 0.0;
    if (xst > 1.0 - tiny) return 1.0;
    return xst;
}


double utl::sqr(const double &number)
{
    return number * number;
}


vecDbl utl::ones(const size_t &n)
{
    vecDbl ones(n);
    for (size_t i = 0u; i < n; ++i)
        ones[i] = 1.0;
    return ones;
}


vecDbl utl::zeros(const size_t &n)
{
    vecDbl zeros(n);
    for (size_t i = 0u; i < n; ++i)
        zeros[i] = 0.0;
    return zeros;
}


vecUns utl::uzeros(const size_t &n)
{
    vecUns zeros(n);
    for (size_t i = 0u; i < n; ++i)
        zeros[i] = 0u;
    return zeros;
}


Matrix utl::matzeros(const size_t &ncol, const size_t &nrow)
{
    Matrix mat(ncol);
    for (size_t i = 0u; i < ncol; ++i)
        mat[i] = zeros(nrow);
    return mat;
}

MatUns utl::matuzeros(const size_t &ncol, const size_t &nrow)
{
    MatUns mat;
    for (size_t i = 0u; i < ncol; ++i)
        mat.push_back(uzeros(nrow));
    return mat;
}

vecDbl utl::rep(const double &x, const size_t &n)
{
    vecDbl reps;
    for (size_t i = 0u; i < n; ++i)
        reps.push_back(x);
    return reps;
}

vecUns utl::repUns(const size_t &x, const size_t &n)
{
    vecUns reps;
    for (size_t i = 0u; i < n; ++i)
        reps.push_back(x);
    return reps;
}


double utl::sum(vecDbl &v)
{
    return std::accumulate(v.begin(), v.end(), 0);
}

double utl::sum(Matrix &m)
{
    double res = 0.0;
    for (auto vec : m)
        res += utl::sum(vec);
    return res;
}

size_t utl::argmin(vecDbl &v)
{
    return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
}

size_t utl::sumu(vecUns &v) {
    size_t sum = 0u;
    for (size_t i = 0u; i < v.size(); ++i)
        sum += v[i];
    return sum;
}

double utl::size2dbl(const size_t &x)
{
    return static_cast<double>(x);
}
