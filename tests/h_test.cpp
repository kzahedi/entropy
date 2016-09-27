#include "h_test.h"

#include <entropy++/Container.h>
#include <entropy++/H.h>
#include <entropy++/sparse/H.h>

#include <iostream>
#include <string>

#include <math.h>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( hTest );


void hTest::testMax()
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

void hTest::testZero()
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
