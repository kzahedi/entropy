#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE pi_test
#include <boost/test/unit_test.hpp>

#include <entropy++/Container.h>
#include <entropy++/PI.h>
#include <entropy++/sparse/PI.h>

#include <iostream>
#include <string>

#include <math.h>

using namespace std;
using namespace entropy;

BOOST_AUTO_TEST_CASE(Sinus)
{

  DContainer container(1000,2);
  for(float i = 0; i < 1000.0; i = i + 1.0)
  {
    container << cos(i/10.0);
    container << sin(i/5.0);
  }

  double **dom = new double*[2];
  dom[0]       = new double[2];
  dom[1]       = new double[2];
  dom[0][0]    = -1.0;
  dom[0][1]    = 1.0;
  dom[1][0]    = -1.0;
  dom[1][1]    = 1.0;

  int *bins    = new int[2];
  bins[0]      = 100;
  bins[1]      = 100;

  container.setDomains(dom);
  container.setBinSizes(bins);

  ULContainer *d  = container.discretise();

  double s = entropy::PI(d);

  BOOST_CHECK_CLOSE(8.196816, s, 0.01); // recalcuate somewhere else

  delete d;
}


BOOST_AUTO_TEST_CASE(SparseVsNonSparse)
{
  DContainer container(1000,2);
  for(float i = 0; i < 1000.0; i = i + 1.0)
  {
    container << cos(i/10.0);
    container << sin(i/5.0);
  }

  double **dom = new double*[2];
  dom[0]       = new double[2];
  dom[1]       = new double[2];
  dom[0][0]    = -1.0;
  dom[0][1]    = 1.0;
  dom[1][0]    = -1.0;
  dom[1][1]    = 1.0;

  int *bins    = new int[2];
  bins[0]      = 100;
  bins[1]      = 100;

  container.setDomains(dom);
  container.setBinSizes(bins);

  ULContainer *d  = container.discretise();

  // PI  pi;
  // PIs pis;

  double s1 = entropy::PI(d);
  double s2 = entropy::sparse::PI(d);
  BOOST_CHECK_CLOSE(s1, s2, 0.01);

  delete d;
}
