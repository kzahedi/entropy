#include "entropy_test.h"

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

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( entropyTest );


void entropyTest::testMax()
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

  CPPUNIT_ASSERT_DOUBLES_EQUAL(log2(1000.0), s, 0.00001);

  double t = entropy::sparse::H(dx);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(log2(1000.0), t, 0.00001);

  delete dx;
}

void entropyTest::testZero()
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

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, s, 0.00001);

  double t = entropy::sparse::H(dx);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, t, 0.00001);

  delete dx;
}

// H(X|Y) = H(X,Y) - H(Y)
void entropyTest::testConditional()
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

  CPPUNIT_ASSERT_DOUBLES_EQUAL(b, a, 0.00001);

  double c = entropy::sparse::ConditionalEntropy(dx, dy);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(a, c, 0.00001);

  delete dx;
  delete dy;
  delete dt;
  delete dz;
}
