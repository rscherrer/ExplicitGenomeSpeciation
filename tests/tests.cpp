#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE testModule

#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(testingABunchOfFunctions)

    BOOST_AUTO_TEST_CASE(OneEqualsOne)
    {
        BOOST_CHECK_EQUAL(1, 1);
    }

BOOST_AUTO_TEST_SUITE_END()