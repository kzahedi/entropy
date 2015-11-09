#include "container_test.h"

#include <entropy/Container.h>

#include <iostream>
#include <string>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( containerTest );


void containerTest::testFilling()
{

  Container container(2,3);

  CPPUNIT_ASSERT_EQUAL(2, container.rows());
  CPPUNIT_ASSERT_EQUAL(3, container.columns());
  
  for(int r = 0; r < container.rows(); r++)
  {
    for(int c = 0; c < container.columns(); c++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, container(r,c), 0.0001);
    }
  }

  for(int r = 0; r < container.rows(); r++)
  {
    for(int c = 0; c < container.columns(); c++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, container.get(r,c), 0.0001);
    }
  }

  for(int i = 0; i < 2 * 3; i++)
  {
    container << i;
  }

  int index = 0;
  for(int r = 0; r < container.rows(); r++)
  {
    for(int c = 0; c < container.columns(); c++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(index, container(r,c), 0.0001);
      index++;
    }
  }

}

void containerTest::testUniformDiscretisationUnary()
{
  Container c(11,1);

  for(int i = 0; i < 11; i++)
  {
    c << ((float)i)/10.0;
  }

  double **domain = new double*[1];
  domain[0]       = new double[2];
  domain[0][0]    = 0.0;
  domain[0][1]    = 1.0;

  int *bins       = new int[1];
  bins[0]         = 10;

  c.setBinSizes(bins);
  c.setDomains(domain);

  Container *d = c.discretise();

  CPPUNIT_ASSERT_EQUAL(0, (int)d->get(0,0));
  CPPUNIT_ASSERT_EQUAL(1, (int)d->get(1,0));
  CPPUNIT_ASSERT_EQUAL(2, (int)d->get(2,0));
  CPPUNIT_ASSERT_EQUAL(3, (int)d->get(3,0));
  CPPUNIT_ASSERT_EQUAL(4, (int)d->get(4,0));
  CPPUNIT_ASSERT_EQUAL(5, (int)d->get(5,0));
  CPPUNIT_ASSERT_EQUAL(6, (int)d->get(6,0));
  CPPUNIT_ASSERT_EQUAL(7, (int)d->get(7,0));
  CPPUNIT_ASSERT_EQUAL(8, (int)d->get(8,0));
  CPPUNIT_ASSERT_EQUAL(9, (int)d->get(9,0));
  CPPUNIT_ASSERT_EQUAL(9, (int)d->get(10,0));

  delete   domain[0];
  delete[] domain;

  delete[] bins;
}

void containerTest::testUniformDiscretisation()
{
  Container c(11,3);

  for(int i = 0; i < 11; i++)
  {
    c << ((float)i)/10.0;
    c << ((float)((i + 1) % 10))/10.0;
    c << ((float)((i + 2) % 10))/10.0;
  }

  double **domain = new double*[3];

  domain[0]       = new double[2];
  domain[0][0]    = 0.0;
  domain[0][1]    = 1.0;

  domain[1]       = new double[2];
  domain[1][0]    = 0.0;
  domain[1][1]    = 1.0;

  domain[2]       = new double[2];
  domain[2][0]    = 0.0;
  domain[2][1]    = 1.0;

  int *bins       = new int[3];
  bins[0]         = 10;
  bins[1]         = 10;
  bins[2]         = 10;

  c.setBinSizes(bins);
  c.setDomains(domain);

  Container *d = c.discretise();

  CPPUNIT_ASSERT_EQUAL(210, (int)d->get(0,0));
  CPPUNIT_ASSERT_EQUAL(321, (int)d->get(1,0));
  CPPUNIT_ASSERT_EQUAL(432, (int)d->get(2,0));
  CPPUNIT_ASSERT_EQUAL(543, (int)d->get(3,0));
  CPPUNIT_ASSERT_EQUAL(654, (int)d->get(4,0));
  CPPUNIT_ASSERT_EQUAL(765, (int)d->get(5,0));
  CPPUNIT_ASSERT_EQUAL(876, (int)d->get(6,0));
  CPPUNIT_ASSERT_EQUAL(987, (int)d->get(7,0));
  CPPUNIT_ASSERT_EQUAL(98,  (int)d->get(8,0));
  CPPUNIT_ASSERT_EQUAL(109, (int)d->get(9,0));
  CPPUNIT_ASSERT_EQUAL(219, (int)d->get(10,0));

  delete domain[0];
  delete[] domain;
}
