#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE matrix_test
#include <boost/test/unit_test.hpp>

#include <entropy++/SparseMatrix.h>
#include <entropy++/Matrix.h>

#include <iostream>
#include <string>
#include <stdlib.h>


using namespace std;


BOOST_AUTO_TEST_CASE(Initialisation)
{
  entropy::SparseMatrix sm;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      BOOST_CHECK_CLOSE(0.0, sm(row,col), 0.001);
    }
  }

  entropy::SparseMatrix sm2(2.0);

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      BOOST_CHECK_CLOSE(2.0, sm2(row,col), 0.001);
    }
  }
}

BOOST_AUTO_TEST_CASE(Set)
{
  entropy::SparseMatrix sm;

  sm(1, 1) = 2.0;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      if(row == 1 && col == 1)
      {
        BOOST_CHECK_CLOSE(2.0, sm(row,col), 0.001);
      }
      else
      {
        BOOST_CHECK_CLOSE(0.0, sm(row,col), 0.001);
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(Add)
{
  time_t t;
  time(&t);
  srand48(t);

  entropy::SparseMatrix sm1;
  entropy::SparseMatrix sm2;

  entropy::Matrix m1(10,10);
  entropy::Matrix m2(10,10);

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
      BOOST_CHECK_CLOSE(m1(row,col), sm1(row,col), 0.001);
    }
  }

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      BOOST_CHECK_CLOSE(m2(row,col), sm2(row,col), 0.001);
    }
  }

  entropy::SparseMatrix sm3 = sm1 + sm2;
  entropy::Matrix        m3 =  m1 +  m2;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      BOOST_CHECK_CLOSE(m3(row,col), sm3(row,col), 0.001);
    }
  }
}

BOOST_AUTO_TEST_CASE(Mul)
{
  time_t t;
  time(&t);
  srand48(t);

  entropy::SparseMatrix sm1;
  entropy::SparseMatrix sm2;

  entropy::Matrix m1(10,10);
  entropy::Matrix m2(10,10);

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
      BOOST_CHECK_CLOSE(m1(row,col), sm1(row,col), 0.001);
    }
  }

  entropy::SparseMatrix sm3 = sm1 * 2.0;
  entropy::Matrix        m3 =  m1 * 2.0;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      BOOST_CHECK_CLOSE(m3(row,col), sm3(row,col), 0.001);
    }
  }
}

BOOST_AUTO_TEST_CASE(Div)
{
  time_t t;
  time(&t);
  srand48(t);

  entropy::SparseMatrix sm1;
  entropy::SparseMatrix sm2;

  entropy::Matrix m1(10,10);
  entropy::Matrix m2(10,10);

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
      BOOST_CHECK_CLOSE(m1(row,col), sm1(row,col), 0.001);
    }
  }

  sm1 /= 2.0;
  m1  /= 2.0;

  for(int row = 0; row < 10; row++)
  {
    for(int col = 0; col < 10; col++)
    {
      BOOST_CHECK_CLOSE(m2(row,col), sm2(row,col), 0.001);
    }
  }
}
