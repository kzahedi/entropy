#include "csv_test.h"

#include <entropy++/Container.h>
#include <entropy++/Csv.h>

#include <fstream>
#include <iostream>
#include <string>

#include <math.h>

using namespace std;
using namespace entropy;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( csvTest );

void csvTest::readTestFile()
{
  Csv *csv = new Csv();
  DContainer *c = csv->read(TEST_CSV, 3, 0, 3, 7);

  CPPUNIT_ASSERT_EQUAL(21, c->rows());
  CPPUNIT_ASSERT_EQUAL(3,  c->columns());

  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, c->get(0,0), 0.00001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0, c->get(0,1), 0.00001);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.0, c->get(0,2), 0.00001);

  for(int i = 1; i < 21; i++)
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL(9+i,  c->get(i,0), 0.00001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(39+i, c->get(i,1), 0.00001);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(79+i, c->get(i,2), 0.00001);
  }
}
