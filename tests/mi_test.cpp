#include "mi_test.h"

#include <entropy++/Container.h>
#include <entropy++/MI.h>
#include <entropy++/MIs.h>

#include <iostream>
#include <string>

#include <math.h>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( miTest );


void miTest::testSinus()
{

  Container X(1000,1);
  Container Y(1000,1);
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

  Container *dx = X.discretise();
  Container *dy = Y.discretise();

  MI mi;
  double s = mi.calculate(dx, dy);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.17221, s, 0.00001); // recalcuate somewhere else

  delete dx;
  delete dy;
}


void miTest::testSparseVsNonSparse()
{
  Container X(1000,1);
  Container Y(1000,1);
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

  Container *dx = X.discretise();
  Container *dy = Y.discretise();

  MI  mi;
  MIs mis;

  double s1 = mi.calculate(dx, dy);
  double s2 = mis.calculate(dx, dy);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(s1, s2, 0.00001);

  delete dx;
  delete dy;
}
