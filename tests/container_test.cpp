#include "container_test.h"

#include <entropy++/Container.h>

#include <iostream>
#include <string>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( containerTest );


void containerTest::testFilling()
{

  DContainer container(2,3);

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

void containerTest::testDropping()
{
  DContainer c(10,3);

  for(int i = 0; i < 30; i++)
  {
    c << i;
  }

  DContainer *d = c.drop(3);

  CPPUNIT_ASSERT_EQUAL(7, d->rows());
  CPPUNIT_ASSERT_EQUAL(3, d->columns());

  int v = 8;
  for(int i = 0; i < d->rows(); i++)
  {
    for(int j = 0; j < d->columns(); j++)
    {
      v++;
      CPPUNIT_ASSERT_EQUAL((int)v, (int)d->get(i,j));
    }
  }

  DContainer *e = d->drop(-3);

  CPPUNIT_ASSERT_EQUAL(4, e->rows());
  CPPUNIT_ASSERT_EQUAL(3, e->columns());

  v = 8;
  for(int i = 0; i < e->rows(); i++)
  {
    for(int j = 0; j < e->columns(); j++)
    {
      v++;
      CPPUNIT_ASSERT_EQUAL((int)v, (int)e->get(i,j));
    }
  }
}

void containerTest::testUniformDiscretisationUnary()
{
  DContainer c(11,1);

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

  ULContainer *d = c.discretise();

  CPPUNIT_ASSERT_EQUAL(0, (int)d->get(0,  0));
  CPPUNIT_ASSERT_EQUAL(1, (int)d->get(1,  0));
  CPPUNIT_ASSERT_EQUAL(2, (int)d->get(2,  0));
  CPPUNIT_ASSERT_EQUAL(3, (int)d->get(3,  0));
  CPPUNIT_ASSERT_EQUAL(4, (int)d->get(4,  0));
  CPPUNIT_ASSERT_EQUAL(5, (int)d->get(5,  0));
  CPPUNIT_ASSERT_EQUAL(6, (int)d->get(6,  0));
  CPPUNIT_ASSERT_EQUAL(7, (int)d->get(7,  0));
  CPPUNIT_ASSERT_EQUAL(8, (int)d->get(8,  0));
  CPPUNIT_ASSERT_EQUAL(9, (int)d->get(9,  0));
  CPPUNIT_ASSERT_EQUAL(9, (int)d->get(10, 0));

  delete   domain[0];
  delete[] domain;
  delete[] bins;
}

void containerTest::testUniformDiscretisation()
{
  DContainer c(11,3);

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

  ULContainer *d = c.discretise();

  CPPUNIT_ASSERT_EQUAL(0,  (int)d->get(0,  0)); // 210
  CPPUNIT_ASSERT_EQUAL(1,  (int)d->get(1,  0)); // 321
  CPPUNIT_ASSERT_EQUAL(2,  (int)d->get(2,  0)); // 432
  CPPUNIT_ASSERT_EQUAL(3,  (int)d->get(3,  0)); // 543
  CPPUNIT_ASSERT_EQUAL(4,  (int)d->get(4,  0)); // 654
  CPPUNIT_ASSERT_EQUAL(5,  (int)d->get(5,  0)); // 765
  CPPUNIT_ASSERT_EQUAL(6,  (int)d->get(6,  0)); // 876
  CPPUNIT_ASSERT_EQUAL(7,  (int)d->get(7,  0)); // 987
  CPPUNIT_ASSERT_EQUAL(8,  (int)d->get(8,  0)); // 98
  CPPUNIT_ASSERT_EQUAL(9,  (int)d->get(9,  0)); // 109
  CPPUNIT_ASSERT_EQUAL(10, (int)d->get(10, 0)); // 210

  delete   domain[0];
  delete   domain[1];
  delete   domain[2];
  delete[] domain;
  delete[] bins;
}

void containerTest::testCopy()
{
  DContainer c(10,3);
  DContainer d(0,0);

  for(int i = 0; i < 30; i++) c << i;

  d = c;

  CPPUNIT_ASSERT_EQUAL(c.rows(),    d.rows());
  CPPUNIT_ASSERT_EQUAL(c.columns(), d.columns());

  for(int i = 0; i < 10; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      CPPUNIT_ASSERT_EQUAL((int)c(i,j), (int)d(i,j));
    }
  }

  c(0, 0) =  10;

  CPPUNIT_ASSERT((int)c(0,0) != (int)d(0,0));
}

void containerTest::testMax()
{
  DContainer c(10,3);

  for(int i = 0; i < 30; i++)
  {
    c << i;
  }

  CPPUNIT_ASSERT_EQUAL(27, (int)c.max(0));
  CPPUNIT_ASSERT_EQUAL(28, (int)c.max(1));
  CPPUNIT_ASSERT_EQUAL(29, (int)c.max(2));
  CPPUNIT_ASSERT_EQUAL(29, (int)c.max());

}

void containerTest::testMin()
{
  DContainer c(10,3);

  for(int i = 0; i < 30; i++)
  {
    c << i;
  }

  CPPUNIT_ASSERT_EQUAL(0, (int)c.min(0));
  CPPUNIT_ASSERT_EQUAL(1, (int)c.min(1));
  CPPUNIT_ASSERT_EQUAL(2, (int)c.min(2));
  CPPUNIT_ASSERT_EQUAL(0, (int)c.min());

}

void containerTest::testExtractColumns()
{
  int rows    = 50;
  int columns = 5;
  DContainer c(rows,columns);
  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < columns; j++)
    {
      int value = (j+1) * 100 + i;
      c << value;
    }
  }

  DContainer *onethreefive = c.columns(3,0,2,4);

  for(int i = 0; i < rows; i++)
  {
    CPPUNIT_ASSERT_EQUAL((int)c(i,0), (int)(*onethreefive)(i,0));
    CPPUNIT_ASSERT_EQUAL((int)c(i,2), (int)(*onethreefive)(i,1));
    CPPUNIT_ASSERT_EQUAL((int)c(i,4), (int)(*onethreefive)(i,2));
  }

}

void containerTest::testNormaliseColumn()
{
  int rows    = 50;
  int columns = 5;
  DContainer c(rows,columns);
  int index = 0;
  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < columns; j++)
    {
      c << index++;
    }
  }

  for(int j = 0; j < columns; j++)
  {
    c.normaliseColumn(j, c.min(j), c.max(j));
  }

  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < columns; j++)
    {
      CPPUNIT_ASSERT(0.0    <= c(i,j));
      CPPUNIT_ASSERT(c(i,j) <= 1.0);
    }
  }
}

void containerTest::testCopyFunc()
{
  int rows    = 50;
  int columns = 5;
  DContainer c(rows,columns);
  int index = 0;
  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < columns; j++)
    {
      c << index++;
    }
  }

  DContainer *d = c.copy();
  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < columns; j++)
    {
      CPPUNIT_ASSERT_EQUAL((int)c(i,j), (int)(*d)(i,j));
    }
  }
}

void containerTest::testMerge()
{
  DContainer c(10,3);
  DContainer d(10,5);

  for(int i = 0; i < 30; i++) c << i;
  for(int i = 0; i < 50; i++) d << i;

  c += d;

  CPPUNIT_ASSERT_EQUAL(10, c.rows());
  CPPUNIT_ASSERT_EQUAL(8,  c.columns());

  int index = 0;
  for(int i = 0; i < 10; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(index++, c(i,j), 0.0001);
    }
  }

  index = 0;
  for(int i = 0; i < 10; i++)
  {
    for(int j = 0; j < 5; j++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(index++, c(i,3+j), 0.0001);
    }
  }
}

void containerTest::testFillMode()
{
  DContainer container(2,3);
  container.setFillMode(FILL_MODE_BY_ROW);

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

  DContainer container2(2,3);
  container2.setFillMode(FILL_MODE_BY_COLUMN);

  for(int i = 0; i < 2 * 3; i++)
  {
    container2 << i;
  }

  index = 0;
  for(int c = 0; c < container2.columns(); c++)
  {
    for(int r = 0; r < container2.rows(); r++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(index, container2(r,c), 0.0001);
      index++;
    }
  }
}


void containerTest::testUnique1()
{
  ULContainer* container = new ULContainer(10,5);

  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;

  ULContainer *unique = container->unique();

  CPPUNIT_ASSERT_EQUAL(1, unique->rows());
}

void containerTest::testUnique2()
{
  ULContainer* container = new ULContainer(10,5);

  *container << 1 << 2 << 3 << 4 << 5;
  *container << 2 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 2 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 1 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 1 << 4 << 5;
  *container << 1 << 2 << 3 << 4 << 5;
  *container << 1 << 2 << 5 << 4 << 5;

  ULContainer *unique = container->unique();

  CPPUNIT_ASSERT_EQUAL(4, unique->rows());
}
