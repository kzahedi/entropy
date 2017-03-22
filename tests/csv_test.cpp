#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE pi_test
#include <boost/test/unit_test.hpp>

#include <entropy++/Container.h>
#include <entropy++/Csv.h>

#include <fstream>
#include <iostream>
#include <string>

#include <math.h>

using namespace std;
using namespace entropy;

BOOST_AUTO_TEST_CASE(readTestFile)
{
  Csv *csv = new Csv();
  DContainer *c = csv->read(TEST_CSV, 3, 0, 3, 7);

  BOOST_CHECK(21 == c->rows());
  BOOST_CHECK(3  == c->columns());

  BOOST_CHECK_CLOSE(0.0, c->get(0,0), 0.00001);
  BOOST_CHECK_CLOSE(3.0, c->get(0,1), 0.00001);
  BOOST_CHECK_CLOSE(7.0, c->get(0,2), 0.00001);

  for(int i = 1; i < 21; i++)
  {
    BOOST_CHECK_CLOSE(9+i,  c->get(i,0), 0.00001);
    BOOST_CHECK_CLOSE(39+i, c->get(i,1), 0.00001);
    BOOST_CHECK_CLOSE(79+i, c->get(i,2), 0.00001);
  }
}
