#include "matrix_test.h"

#include <entropy++/SparseMatrix.h>
#include <entropy++/Matrix.h>

#include <iostream>
#include <string>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( matrixTest );

void matrixTest::testInitialisation()
{
  SparseMatrix sm;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, sm(row,col), 0.0001);
    }
  }

  SparseMatrix sm2(2.0);

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, sm2(row,col), 0.0001);
    }
  }
}

void matrixTest::testSet()
{
  SparseMatrix sm;

  sm(1, 1) = 2.0;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      if(row == 1 && col == 1)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0, sm(row,col), 0.0001);
      }
      else
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, sm(row,col), 0.0001);
      }
    }
  }
}

void matrixTest::testAdd()
{
  time_t t;
  time(&t);
  srand48(t);

  SparseMatrix sm1;
  SparseMatrix sm2;

  Matrix m1(10,10);
  Matrix m2(10,10);

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      double r1 = drand48();
      double r2 = drand48();
      m1(row, col)  = r1;
      sm1(row, col) = r1;
      m2(row, col)  = r2;
      sm2(row, col) = r2;
    }
  }

  m1  += m2;
  sm1 += sm2;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m1(row,col), sm1(row,col), 0.0001);
    }
  }

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m2(row,col), sm2(row,col), 0.0001);
    }
  }

  SparseMatrix sm3 = sm1 + sm2;
  Matrix        m3 =  m1 +  m2;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m3(row,col), sm3(row,col), 0.0001);
    }
  }
}

void matrixTest::testMul()
{
  time_t t;
  time(&t);
  srand48(t);

  SparseMatrix sm1;
  SparseMatrix sm2;

  Matrix m1(10,10);
  Matrix m2(10,10);

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      double r1 = drand48();
      double r2 = drand48();
      m1(row, col)  = r1;
      sm1(row, col) = r1;
      m2(row, col)  = r2;
      sm2(row, col) = r2;
    }
  }

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m1(row,col), sm1(row,col), 0.0001);
    }
  }

  SparseMatrix sm3 = sm1 * 2.0;
  Matrix        m3 =  m1 * 2.0;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m3(row,col), sm3(row,col), 0.0001);
    }
  }
}

void matrixTest::testDiv()
{
  time_t t;
  time(&t);
  srand48(t);

  SparseMatrix sm1;
  SparseMatrix sm2;

  Matrix m1(10,10);
  Matrix m2(10,10);

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      double r1 = drand48();
      double r2 = drand48();
      m1(row, col)  = r1;
      sm1(row, col) = r1;
      m2(row, col)  = r2;
      sm2(row, col) = r2;
    }
  }

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m1(row,col), sm1(row,col), 0.0001);
    }
  }

  sm1 /= 2.0;
  m1  /= 2.0;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(m2(row,col), sm2(row,col), 0.0001);
    }
  }
}
