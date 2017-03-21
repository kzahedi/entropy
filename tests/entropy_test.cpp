#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE entropy_test
#include <boost/test/unit_test.hpp>

#include <entropy++/Container.h>

#include <entropy++/H.h>
#include <entropy++/sparse/H.h>

#include <entropy++/ConditionalEntropy.h>
#include <entropy++/sparse/ConditionalEntropy.h>

#include <iostream>
#include <string>
#include <stdlib.h>

#include <math.h>

using namespace std;
using namespace entropy;


BOOST_AUTO_TEST_CASE(Max)
{
  DContainer X(1000,1);
  for(float i = 0; i < 1000.0; i = i + 1.0)
  {
    X << i;
  }

  double **dom = new double*[1];
  dom[0]       = new double[2];
  dom[0][0]    = 0.0;
  dom[0][1]    = 999.0;

  int *bins    = new int[1];
  bins[0]      = 1000;

  X.setDomains(dom);
  X.setBinSizes(bins);

  ULContainer *dx = X.discretise();

  double s = entropy::H(dx);

  BOOST_CHECK_CLOSE(log2(1000.0), s, 0.001);

  double t = entropy::sparse::H(dx);

  BOOST_CHECK_CLOSE(log2(1000.0), t, 0.001);

  delete dx;
}

BOOST_AUTO_TEST_CASE(Zero)
{
  DContainer X(1000,1);
  for(float i = 0; i < 1000.0; i = i + 1.0)
  {
    X << 0.0;
  }

  double **dom = new double*[1];
  dom[0]       = new double[2];
  dom[0][0]    = 0.0;
  dom[0][1]    = 999.0;

  int *bins    = new int[1];
  bins[0]      = 1000;

  X.setDomains(dom);
  X.setBinSizes(bins);

  ULContainer *dx = X.discretise();

  double s = entropy::H(dx);

  BOOST_CHECK_CLOSE(0.0, s, 0.001);

  double t = entropy::sparse::H(dx);

  BOOST_CHECK_CLOSE(0.0, t, 0.001);

  delete dx;
}

// H(X|Y) = H(X,Y) - H(Y)
BOOST_AUTO_TEST_CASE(Conditional)
{
  DContainer X(1000,1);
  DContainer Y(1000,1);

#ifdef __APPLE__
  sranddev();
#else // __APPLE__
  srand(time(NULL));
#endif // __APPLE__

  for(int i = 0; i < 1000; i++)
  {
    X << ((float)rand())/((float)RAND_MAX);
    Y << ((float)rand())/((float)RAND_MAX);
  }

  double **dom = new double*[1];
  dom[0]       = new double[2];
  dom[0][0]    = 0.0;
  dom[0][1]    = 1.0;

  int *bins    = new int[1];
  bins[0]      = 1000;

  X.setDomains(dom);
  X.setBinSizes(bins);

  Y.setDomains(dom);
  Y.setBinSizes(bins);

  ULContainer *dx = X.discretise();
  ULContainer *dy = Y.discretise();
  ULContainer *dt = dx->copy();
  *dt += *dy;
  ULContainer *dz = dt->combineDiscretisedColumns();

  double a = entropy::ConditionalEntropy(dx, dy);

  double b = entropy::H(dz) - entropy::H(dy);

  BOOST_CHECK_CLOSE(b, a, 0.001);

  double c = entropy::sparse::ConditionalEntropy(dx, dy);

  BOOST_CHECK_CLOSE(a, c, 0.001);

  delete dx;
  delete dy;
  delete dt;
  delete dz;
}
