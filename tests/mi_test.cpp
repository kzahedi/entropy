#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mi_test
#include <boost/test/unit_test.hpp>

#include <entropy++/Container.h>
#include <entropy++/MI.h>
#include <entropy++/sparse/MI.h>

#include <iostream>
#include <string>

#include <math.h>

using namespace std;
using namespace entropy;


BOOST_AUTO_TEST_CASE(Sinus)
{

  DContainer X(1000,1);
  DContainer Y(1000,1);
  for(float i = 0; i < 1000.0; i = i + 1.0)
  {
    X << cos(i/10.0);
    Y << sin(i/5.0);
  }

  double **dom = new double*[1];
  dom[0]       = new double[2];
  dom[0][0]    = -1.0;
  dom[0][1]    = 1.0;

  int *bins    = new int[1];
  bins[0]      = 100;

  X.setDomains(dom);
  X.setBinSizes(bins);

  Y.setDomains(dom);
  Y.setBinSizes(bins);

  ULContainer *dx = X.discretise();
  ULContainer *dy = Y.discretise();

  double s = entropy::MI(dx, dy);

  BOOST_CHECK_CLOSE(4.17221, s, 0.001); // recalcuate somewhere else

  delete dx;
  delete dy;
}


BOOST_AUTO_TEST_CASE(SparseVsNonSparse)
{
  DContainer X(1000,1);
  DContainer Y(1000,1);
  for(float i = 0; i < 1000.0; i = i + 1.0)
  {
    X << cos(i/10.0);
    Y << sin(i/5.0);
  }

  double **dom = new double*[1];
  dom[0]       = new double[2];
  dom[0][0]    = -1.0;
  dom[0][1]    = 1.0;

  int *bins    = new int[1];
  bins[0]      = 100;

  X.setDomains(dom);
  X.setBinSizes(bins);

  Y.setDomains(dom);
  Y.setBinSizes(bins);

  ULContainer *dx = X.discretise();
  ULContainer *dy = Y.discretise();

  double s1 = entropy::MI(dx, dy);
  double s2 = entropy::sparse::MI(dx, dy);
  BOOST_CHECK_CLOSE(s1, s2, 0.001);

  delete dx;
  delete dy;
}
