#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "utils.h"

BOOST_AUTO_TEST_CASE(testSqr)
{
    BOOST_CHECK_EQUAL(sqr(2.0), 4.0);
    BOOST_CHECK_EQUAL(sqr(4.0), 16.0);
}

/*
BOOST_AUTO_TEST_CASE(test_calc_mean)
{
    const double measured{
            calc_mean( {1.0, 2.0, 3.0} )
    };
    const double expected{2.0};
    BOOST_CHECK_EQUAL(measured, expected);
}

BOOST_AUTO_TEST_CASE(test_calc_mean_needs_nonempty_vector)
{
    std::vector<double> empty;
    BOOST_CHECK_THROW(
            calc_mean(empty),
            std::invalid_argument
    );
}

 */