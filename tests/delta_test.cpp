#include "delta_test.h"

#include <entropy++/iterativescaling/Delta.h>

#include <iostream>
#include <vector>

#include <math.h>

using namespace std;
using namespace entropy::iterativescaling;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( deltaTest );


void deltaTest::testDeltaMatch()
{
  vector<unsigned long> xrow;
  xrow.push_back(0);
  xrow.push_back(1);
  xrow.push_back(2);
  xrow.push_back(3);
  xrow.push_back(4);
  xrow.push_back(5);
  xrow.push_back(6);
  xrow.push_back(7);
  xrow.push_back(8);
  xrow.push_back(9);

  vector<unsigned long> xrowmatch;
  xrowmatch.push_back(2);
  xrowmatch.push_back(4);
  xrowmatch.push_back(6);
  xrowmatch.push_back(8);

  vector<int> xIndices;
  xIndices.push_back(2);
  xIndices.push_back(4);
  xIndices.push_back(6);
  xIndices.push_back(8);

  vector<unsigned long> yrow;
  yrow.push_back(0);
  yrow.push_back(1);
  yrow.push_back(2);
  yrow.push_back(3);
  yrow.push_back(4);
  yrow.push_back(5);
  yrow.push_back(6);
  yrow.push_back(7);
  yrow.push_back(8);
  yrow.push_back(9);

  vector<unsigned long> yrowmatch;
  yrowmatch.push_back(1);
  yrowmatch.push_back(3);
  yrowmatch.push_back(5);
  yrowmatch.push_back(7);

  vector<int> yIndices;
  yIndices.push_back(1);
  yIndices.push_back(3);
  yIndices.push_back(5);
  yIndices.push_back(7);

  Delta *d = new Delta(xrow, xIndices, yrow, yIndices);
  CPPUNIT_ASSERT_EQUAL(true, d->match(xrowmatch, yrowmatch));
}

void deltaTest::testDeltaMatchXY()
{
  vector<unsigned long> xrow;
  xrow.push_back(0);
  xrow.push_back(1);
  xrow.push_back(2);
  xrow.push_back(3);
  xrow.push_back(4);
  xrow.push_back(5);
  xrow.push_back(6);
  xrow.push_back(7);
  xrow.push_back(8);
  xrow.push_back(9);

  vector<unsigned long> xrow2;
  xrow2.push_back(1);
  xrow2.push_back(2);
  xrow2.push_back(3);
  xrow2.push_back(4);
  xrow2.push_back(5);
  xrow2.push_back(6);
  xrow2.push_back(7);
  xrow2.push_back(8);
  xrow2.push_back(9);
  xrow2.push_back(0);

  vector<unsigned long> xrowmatch;
  xrowmatch.push_back(0);
  xrowmatch.push_back(1);
  xrowmatch.push_back(2);
  xrowmatch.push_back(3);
  xrowmatch.push_back(4);
  xrowmatch.push_back(5);
  xrowmatch.push_back(6);
  xrowmatch.push_back(7);
  xrowmatch.push_back(8);
  xrowmatch.push_back(9);

  vector<int> xIndices;
  xIndices.push_back(2);
  xIndices.push_back(4);
  xIndices.push_back(6);
  xIndices.push_back(8);

  vector<unsigned long> yrow;
  yrow.push_back(0);
  yrow.push_back(1);
  yrow.push_back(2);
  yrow.push_back(3);
  yrow.push_back(4);
  yrow.push_back(5);
  yrow.push_back(6);
  yrow.push_back(7);
  yrow.push_back(8);
  yrow.push_back(9);

  vector<unsigned long> yrow2;
  yrow2.push_back(1);
  yrow2.push_back(2);
  yrow2.push_back(3);
  yrow2.push_back(4);
  yrow2.push_back(5);
  yrow2.push_back(6);
  yrow2.push_back(7);
  yrow2.push_back(8);
  yrow2.push_back(9);
  yrow2.push_back(0);

  vector<unsigned long> yrowmatch;
  yrowmatch.push_back(0);
  yrowmatch.push_back(1);
  yrowmatch.push_back(2);
  yrowmatch.push_back(3);
  yrowmatch.push_back(4);
  yrowmatch.push_back(5);
  yrowmatch.push_back(6);
  yrowmatch.push_back(7);
  yrowmatch.push_back(8);
  yrowmatch.push_back(9);

  vector<int> yIndices;
  yIndices.push_back(1);
  yIndices.push_back(3);
  yIndices.push_back(5);
  yIndices.push_back(7);

  Delta *d1 = new Delta(xrow,  xIndices, yrow,  yIndices);
  Delta *d2 = new Delta(xrow2, xIndices, yrow2, yIndices);
  CPPUNIT_ASSERT_EQUAL(true,  d1->matchXY(xrowmatch, yrowmatch));
  CPPUNIT_ASSERT_EQUAL(false, d2->matchXY(xrowmatch, yrowmatch));
}
