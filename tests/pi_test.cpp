#include "pi_test.h"

#include <entropy++/Container.h>
#include <entropy++/PI.h>
#include <entropy++/PIs.h>

#include <iostream>
#include <string>

#include <math.h>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( piTest );


void piTest::testSinus()
{

  Container container(1000,2);
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

  Container *d  = container.discretise();

  PI pi;

  double s = pi.calculate(d);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(9.56171, s, 0.00001); // recalcuate somewhere else

  delete d;
}


void piTest::testSparseVsNonSparse()
{
  Container container(1000,2);
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

  Container *d  = container.discretise();

  PI  pi;
  PIs pis;

  double s1 = pi.calculate(d);
  double s2 = pis.calculate(d);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(s1, s2, 0.00001);

  delete d;
}
