#define BOOST_TEST_MODULE omp_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <stdio.h>
#include <iostream>
#include <omp.h>

using namespace std;

BOOST_AUTO_TEST_SUITE(OpenMP)

BOOST_AUTO_TEST_CASE(OMP)
{
  vector<int> s(10);
#pragma omp parallel for
  for(int i = 0; i < 100; i++)
  {
    printf("Hello from thread %d, nthreads %d -> %d\n", omp_get_thread_num(), omp_get_num_threads(),i);
    s[i % 10] += i;
  }

  for(vector<int>::iterator i = s.begin(); i != s.end(); i++)
  {
    cout << *i << endl;
  }
}
BOOST_AUTO_TEST_SUITE_END()
