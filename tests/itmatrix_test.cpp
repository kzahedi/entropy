#include "itmatrix_test.h"

#include <entropy++/Container.h>
#include <entropy++/MI.h>
#include <entropy++/sparse/MI.h>

#include <entropy++/iterativescaling/ITMatrix.h>

#include <iostream>
#include <string>

using namespace std;
using namespace entropy;
using namespace entropy::iterativescaling;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(itMatrixTest);

void itMatrixTest::testFillXBinary()
{

  ULContainer *xAlphabet = new ULContainer(2,1); // alphabet
  *xAlphabet << 0 << 1;
  ULContainer *yAlphabet = new ULContainer(2,1); // alphabet
  *yAlphabet << 0 << 1;

  ivvector alphX(1,ivector(0));
  alphX[0].push_back(0);

  ivvector alphY(1,ivector(0));
  alphY[0].push_back(0);

  dvector lambda(1);
  lambda[0] = 1;

  int xColumns      = 2;
  int xRows         = 1000;
  int yColumns      = 1;

  ULContainer* xData = new ULContainer(xRows,xColumns);

  for(int i = 0;i < xRows; i++)
  {
    for(int j = 0;j < xColumns; j++)
    {
      double z = rand() % xAlphabet->rows(); // random indices for the alphabet
      *xData << (*xAlphabet)(z,0);
    }
  }

  ULContainer* yData = new ULContainer(xRows,yColumns);

  for(int i = 0;i < xRows; i++)
  {
    for(int j = 0;j < yColumns; j++)
    {
      double z = rand() % yAlphabet->rows(); // random indices for the alphabet
      *yData << (*yAlphabet)(z,0);
    }
  }


  ITMatrix *it = new ITMatrix(xData,
                              yData,
                              xAlphabet,
                              yAlphabet,
                              alphX,
                              alphY,
                              0.0);

  it->fillX();
  // 0 0
  // 1 0
  // 0 1
  // 1 1
  CPPUNIT_ASSERT_EQUAL(0, it->getFillX(0,0));
  CPPUNIT_ASSERT_EQUAL(0, it->getFillX(0,1));

  CPPUNIT_ASSERT_EQUAL(1, it->getFillX(1,0));
  CPPUNIT_ASSERT_EQUAL(0, it->getFillX(1,1));

  CPPUNIT_ASSERT_EQUAL(0, it->getFillX(2,0));
  CPPUNIT_ASSERT_EQUAL(1, it->getFillX(2,1));

  CPPUNIT_ASSERT_EQUAL(1, it->getFillX(3,0));
  CPPUNIT_ASSERT_EQUAL(1, it->getFillX(3,1));
}

void itMatrixTest::testFillXMore()
{

  ULContainer *xAlphabet = new ULContainer(4,1); // alphabet
  *xAlphabet << 0 << 1 << 2 << 3;
  ULContainer *yAlphabet = new ULContainer(2,1); // alphabet
  *yAlphabet << 0 << 1;

  ivvector alphX(1,ivector(0));
  alphX[0].push_back(0);

  ivvector alphY(1,ivector(0));
  alphY[0].push_back(0);

  dvector lambda(1);
  lambda[0] = 1;

  int xColumns      = 2;
  int xRows         = 1000;
  int yColumns      = 1;

  ULContainer* xData = new ULContainer(xRows,xColumns);

  for(int i = 0;i < xRows; i++)
  {
    for(int j = 0;j < xColumns; j++)
    {
      double z = rand() % xAlphabet->rows(); // random indices for the alphabet
      *xData << (*xAlphabet)(z,0);
    }
  }

  ULContainer* yData = new ULContainer(xRows,yColumns);

  for(int i = 0;i < xRows; i++)
  {
    for(int j = 0;j < yColumns; j++)
    {
      double z = rand() % yAlphabet->rows(); // random indices for the alphabet
      *yData << (*yAlphabet)(z,0);
    }
  }


  ITMatrix *it = new ITMatrix(xData,
                              yData,
                              xAlphabet,
                              yAlphabet,
                              alphX,
                              alphY,
                              0.0);

  it->fillX();
  // 0 0
  // 1 0
  // 2 0
  // 3 0

  // 0 1
  // 1 1
  // 2 1
  // 3 1

  // 0 2
  // 1 2
  // 2 2
  // 3 2

  // 0 3
  // 1 3
  // 2 3
  // 3 3

  for(int i = 0; i < 4; i++)
  {
    CPPUNIT_ASSERT_EQUAL(0, it->getFillX(i*4+0,0));
    CPPUNIT_ASSERT_EQUAL(1, it->getFillX(i*4+1,0));
    CPPUNIT_ASSERT_EQUAL(2, it->getFillX(i*4+2,0));
    CPPUNIT_ASSERT_EQUAL(3, it->getFillX(i*4+3,0));

    CPPUNIT_ASSERT_EQUAL(i, it->getFillX(i*4+0,1));
    CPPUNIT_ASSERT_EQUAL(i, it->getFillX(i*4+1,1));
    CPPUNIT_ASSERT_EQUAL(i, it->getFillX(i*4+2,1));
    CPPUNIT_ASSERT_EQUAL(i, it->getFillX(i*4+3,1));
  }
}
