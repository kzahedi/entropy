#define BOOST_TEST_MODULE omp_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>

#include <omp.h>

using namespace std;

BOOST_AUTO_TEST_SUITE(OpenMP)

BOOST_AUTO_TEST_CASE(OMP)
{
#pragma omp parallel
  {
    int threads = omp_get_num_threads();
    int id = omp_get_thread_num();
    cout << "hello from thread: " << id << " out of " << threads << endl;
    sleep(2);
  }
}
BOOST_AUTO_TEST_SUITE_END()
