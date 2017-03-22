#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE container_test
#include <boost/test/unit_test.hpp>

#include <entropy++/Container.h>

#include <iostream>
#include <string>

using namespace std;
using namespace entropy;

BOOST_AUTO_TEST_SUITE(First)

BOOST_AUTO_TEST_CASE(Filling)
{
  DContainer container(2,3);

  BOOST_CHECK(2 == container.rows());
  BOOST_CHECK(3 == container.columns());

  for(int r = 0; r < container.rows(); r++)
  {
    for(int c = 0; c < container.columns(); c++)
    {
      BOOST_CHECK_CLOSE(0.0, container(r,c), 0.001);
    }
  }

  for(int r = 0; r < container.rows(); r++)
  {
    for(int c = 0; c < container.columns(); c++)
    {
      BOOST_CHECK_CLOSE(0.0, container.get(r,c), 0.001);
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
      BOOST_CHECK_CLOSE(index, container(r,c), 0.001);
      index++;
    }
  }

}

BOOST_AUTO_TEST_CASE(Dropping)
{
  DContainer c(10,3);

  for(int i = 0; i < 30; i++)
  {
    c << i;
  }

  DContainer *d = c.drop(3);

  BOOST_CHECK(7 == d->rows());
  BOOST_CHECK(3 == d->columns());

  int v = 8;
  for(int i = 0; i < d->rows(); i++)
  {
    for(int j = 0; j < d->columns(); j++)
    {
      v++;
      BOOST_CHECK((int)v == (int)d->get(i,j));
    }
  }

  DContainer *e = d->drop(-3);

  BOOST_CHECK(4 == e->rows());
  BOOST_CHECK(3 == e->columns());

  v = 8;
  for(int i = 0; i < e->rows(); i++)
  {
    for(int j = 0; j < e->columns(); j++)
    {
      v++;
      BOOST_CHECK((int)v == (int)e->get(i,j));
    }
  }
}

BOOST_AUTO_TEST_CASE(UniformDiscretisationUnary)
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

  BOOST_CHECK(0 == (int)d->get(0,  0));
  BOOST_CHECK(1 == (int)d->get(1,  0));
  BOOST_CHECK(2 == (int)d->get(2,  0));
  BOOST_CHECK(3 == (int)d->get(3,  0));
  BOOST_CHECK(4 == (int)d->get(4,  0));
  BOOST_CHECK(5 == (int)d->get(5,  0));
  BOOST_CHECK(6 == (int)d->get(6,  0));
  BOOST_CHECK(7 == (int)d->get(7,  0));
  BOOST_CHECK(8 == (int)d->get(8,  0));
  BOOST_CHECK(9 == (int)d->get(9,  0));
  BOOST_CHECK(9 == (int)d->get(10, 0));

  delete   domain[0];
  delete[] domain;
  delete[] bins;
}

BOOST_AUTO_TEST_CASE(UniformDiscretisationByColumn)
{
  int n = 10;
  DContainer c(n,2);

  for(int i = 0; i < n; i++)
  {
    double j = ((double)i)/((float)(n-1));
    double k = ((float)((i+5)%n))/((float)(n-1));
    c << 2.0 * j - 1.0;
    c << 2.0 * k - 1.0;
  }


  double **domain = new double*[1];
  domain[0]       = new double[2];
  domain[0][0]    = -1.0;
  domain[0][1]    =  1.0;

  int *bins       = new int[1];
  bins[0]         = 10;

  c.setBinSizes(10);
  c.setDomains(-1.0, 1.0);

  ULContainer *d = c.discretiseByColumn(false);

  BOOST_CHECK(0 == (int)d->get(0, 0));
  BOOST_CHECK(1 == (int)d->get(1, 0));
  BOOST_CHECK(2 == (int)d->get(2, 0));
  BOOST_CHECK(3 == (int)d->get(3, 0));
  BOOST_CHECK(4 == (int)d->get(4, 0));
  BOOST_CHECK(5 == (int)d->get(5, 0));
  BOOST_CHECK(6 == (int)d->get(6, 0));
  BOOST_CHECK(7 == (int)d->get(7, 0));
  BOOST_CHECK(8 == (int)d->get(8, 0));
  BOOST_CHECK(9 == (int)d->get(9, 0));

  BOOST_CHECK(5 == (int)d->get(0, 1));
  BOOST_CHECK(6 == (int)d->get(1, 1));
  BOOST_CHECK(7 == (int)d->get(2, 1));
  BOOST_CHECK(8 == (int)d->get(3, 1));
  BOOST_CHECK(9 == (int)d->get(4, 1));
  BOOST_CHECK(0 == (int)d->get(5, 1));
  BOOST_CHECK(1 == (int)d->get(6, 1));
  BOOST_CHECK(2 == (int)d->get(7, 1));
  BOOST_CHECK(3 == (int)d->get(8, 1));
  BOOST_CHECK(4 == (int)d->get(9, 1));
}

BOOST_AUTO_TEST_CASE(UniformDiscretisationUnary2)
{
  int n = 10;
  DContainer c(n+1,1);

  for(int i = 0; i < n+1; i++)
  {
    c << 2.0 * ((float)i)/((float)n) - 1.0;
  }

  double **domain = new double*[1];
  domain[0]       = new double[2];
  domain[0][0]    = -1.0;
  domain[0][1]    =  1.0;

  int *bins       = new int[1];
  bins[0]         = 20;

  c.setBinSizes(bins);
  c.setDomains(domain);

  ULContainer *d = c.discretise();

  BOOST_CHECK(0 == (int)d->get(0,  0));
  BOOST_CHECK(1 == (int)d->get(1,  0));
  BOOST_CHECK(2 == (int)d->get(2,  0));
  BOOST_CHECK(3 == (int)d->get(3,  0));
  BOOST_CHECK(4 == (int)d->get(4,  0));
  BOOST_CHECK(5 == (int)d->get(5,  0));
  BOOST_CHECK(6 == (int)d->get(6,  0));
  BOOST_CHECK(7 == (int)d->get(7,  0));
  BOOST_CHECK(8 == (int)d->get(8,  0));
  BOOST_CHECK(9 == (int)d->get(9,  0));

  delete   domain[0];
  delete[] domain;
  delete[] bins;
}

BOOST_AUTO_TEST_CASE(UniformDiscretisation)
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

  BOOST_CHECK(0  == (int)d->get(0,  0)); // 210
  BOOST_CHECK(1  == (int)d->get(1,  0)); // 321
  BOOST_CHECK(2  == (int)d->get(2,  0)); // 432
  BOOST_CHECK(3  == (int)d->get(3,  0)); // 543
  BOOST_CHECK(4  == (int)d->get(4,  0)); // 654
  BOOST_CHECK(5  == (int)d->get(5,  0)); // 765
  BOOST_CHECK(6  == (int)d->get(6,  0)); // 876
  BOOST_CHECK(7  == (int)d->get(7,  0)); // 987
  BOOST_CHECK(8  == (int)d->get(8,  0)); // 98
  BOOST_CHECK(9  == (int)d->get(9,  0)); // 109
  BOOST_CHECK(10 == (int)d->get(10, 0)); // 210

  delete   domain[0];
  delete   domain[1];
  delete   domain[2];
  delete[] domain;
  delete[] bins;
}

BOOST_AUTO_TEST_CASE(Copy)
{
  DContainer c(10,3);
  DContainer d(0,0);

  for(int i = 0; i < 30; i++) c << i;

  d = c;

  BOOST_CHECK(c.rows()    == d.rows());
  BOOST_CHECK(c.columns() == d.columns());

  for(int i = 0; i < 10; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      BOOST_CHECK((int)c(i,j) == (int)d(i,j));
    }
  }

  c(0, 0) =  10;

  BOOST_CHECK((int)c(0,0) != (int)d(0,0));
}

BOOST_AUTO_TEST_CASE(Max)
{
  DContainer c(10,3);

  for(int i = 0; i < 30; i++)
  {
    c << i;
  }

  BOOST_CHECK(27 == (int)c.max(0));
  BOOST_CHECK(28 == (int)c.max(1));
  BOOST_CHECK(29 == (int)c.max(2));
  BOOST_CHECK(29 == (int)c.max());

}

BOOST_AUTO_TEST_CASE(Min)
{
  DContainer c(10,3);

  for(int i = 0; i < 30; i++) c << i;

  BOOST_CHECK(0 == (int)c.min(0));
  BOOST_CHECK(1 == (int)c.min(1));
  BOOST_CHECK(2 == (int)c.min(2));
  BOOST_CHECK(0 == (int)c.min());

}

BOOST_AUTO_TEST_CASE(ExtractColumns)
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
    BOOST_CHECK((int)c(i,0) == (int)(*onethreefive)(i,0));
    BOOST_CHECK((int)c(i,2) == (int)(*onethreefive)(i,1));
    BOOST_CHECK((int)c(i,4) == (int)(*onethreefive)(i,2));
  }

}

BOOST_AUTO_TEST_CASE(NormaliseColumn)
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
      BOOST_CHECK((0.0    <= c(i,j)));
      BOOST_CHECK((c(i,j) <= 1.0));
    }
  }
}

BOOST_AUTO_TEST_CASE(CopyFunc)
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
      BOOST_CHECK((int)c(i,j) == (int)(*d)(i,j));
    }
  }
}

BOOST_AUTO_TEST_CASE(Merge)
{
  DContainer c(10,3);
  DContainer d(10,5);

  for(int i = 0; i < 30; i++) c << i;
  for(int i = 0; i < 50; i++) d << i;

  c += d;

  BOOST_CHECK(10 == c.rows());
  BOOST_CHECK(8  == c.columns());

  int index = 0;
  for(int i = 0; i < 10; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      BOOST_CHECK_CLOSE(index++, c(i,j), 0.001);
    }
  }

  index = 0;
  for(int i = 0; i < 10; i++)
  {
    for(int j = 0; j < 5; j++)
    {
      BOOST_CHECK_CLOSE(index++, c(i,3+j), 0.001);
    }
  }
}

BOOST_AUTO_TEST_CASE(FillMode)
{
  DContainer container(2,3);
  container.setFillMode(FILL_MODE_BY_ROW);

  BOOST_CHECK(2 == container.rows());
  BOOST_CHECK(3 == container.columns());

  for(int r = 0; r < container.rows(); r++)
  {
    for(int c = 0; c < container.columns(); c++)
    {
      BOOST_CHECK_CLOSE(0.0, container(r,c), 0.001);
    }
  }

  for(int r = 0; r < container.rows(); r++)
  {
    for(int c = 0; c < container.columns(); c++)
    {
      BOOST_CHECK_CLOSE(0.0, container.get(r,c), 0.001);
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
      BOOST_CHECK_CLOSE(index, container(r,c), 0.001);
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
      BOOST_CHECK_CLOSE(index, container2(r,c), 0.001);
      index++;
    }
  }
}


BOOST_AUTO_TEST_CASE(Unique1)
{
  IContainer* container  = new IContainer(10,5);
  IContainer* uniqueTest = new IContainer(1,5);

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

  *uniqueTest << 1 << 2 << 3 << 4 << 5;

  IContainer *unique = container->unique();

  BOOST_CHECK(1 == unique->rows());
  BOOST_CHECK(uniqueTest->equals(unique));
}

BOOST_AUTO_TEST_CASE(Unique2)
{
  IContainer* container  = new IContainer(10,5);
  IContainer* uniqueTest = new IContainer(4,5);

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

  *uniqueTest << 1 << 2 << 3 << 4 << 5;
  *uniqueTest << 2 << 2 << 3 << 4 << 5;
  *uniqueTest << 1 << 2 << 1 << 4 << 5;
  *uniqueTest << 1 << 2 << 5 << 4 << 5;

  IContainer *unique = container->unique();

  BOOST_CHECK(uniqueTest->equals(unique));
  BOOST_CHECK(4 == unique->rows());
}

BOOST_AUTO_TEST_CASE(Find1)
{
  IContainer* container = new IContainer(10,5);

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

  int *v = new int[5];
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;
  v[3] = 4;
  v[4] = 5;
  int r = container->find(v);

  BOOST_CHECK(0 == r);
}

BOOST_AUTO_TEST_CASE(Find2)
{
  IContainer* container = new IContainer(10,5);

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

  int *v = new int[5];
  int r = -1;

  v[0] = 1; v[1] = 2; v[2] = 3; v[3] = 4; v[4] = 5;
  r = container->find(v);
  BOOST_CHECK(0 == r);

  v[0] = 2; v[1] = 2; v[2] = 3; v[3] = 4; v[4] = 5;
  r = container->find(v);
  BOOST_CHECK(1 == r);

  v[0] = 1; v[1] = 2; v[2] = 1; v[3] = 4; v[4] = 5;
  r = container->find(v);
  BOOST_CHECK(5 == r);

  v[0] = 1; v[1] = 2; v[2] = 5; v[3] = 4; v[4] = 5;
  r = container->find(v);
  BOOST_CHECK(9 == r);
}

BOOST_AUTO_TEST_CASE(FindList1)
{
  IContainer* container = new IContainer(10,5);

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

  int *v = new int[5];
  v[0] = 1;
  v[1] = 2;
  v[2] = 3;
  v[3] = 4;
  v[4] = 5;
  vector<int> r = container->findlist(v);

  BOOST_CHECK(10 == (int)r.size());
  for(int i = 0; i < 10; i++)
  {
    BOOST_CHECK(i == r[i]);
  }
}

BOOST_AUTO_TEST_CASE(FindList2)
{
  IContainer* container = new IContainer(10,5);

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

  int *v = new int[5];
  vector<int> r;

  v[0] = 1; v[1] = 2; v[2] = 3; v[3] = 4; v[4] = 5;
  r = container->findlist(v);
  BOOST_CHECK(5 == (int)r.size());
  BOOST_CHECK(0 == r[0]);
  BOOST_CHECK(2 == r[1]);
  BOOST_CHECK(4 == r[2]);
  BOOST_CHECK(6 == r[3]);
  BOOST_CHECK(8 == r[4]);

  v[0] = 2; v[1] = 2; v[2] = 3; v[3] = 4; v[4] = 5;
  r = container->findlist(v);
  BOOST_CHECK(2 == (int)r.size());
  BOOST_CHECK(1 == r[0]);
  BOOST_CHECK(3 == r[1]);

  v[0] = 1; v[1] = 2; v[2] = 1; v[3] = 4; v[4] = 5;
  r = container->findlist(v);
  BOOST_CHECK(2 == (int)r.size());
  BOOST_CHECK(5 == r[0]);
  BOOST_CHECK(7 == r[1]);

  v[0] = 1; v[1] = 2; v[2] = 5; v[3] = 4; v[4] = 5;
  r = container->findlist(v);
  BOOST_CHECK(1 == (int)r.size());
  BOOST_CHECK(9 == r[0]);
}


BOOST_AUTO_TEST_CASE(Find1ByContainer)
{
  IContainer* container = new IContainer(10,5);

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

  IContainer* unique = new IContainer(1,5);
  *unique << 1 << 2 << 3 << 4 << 5;

  int r = container->find(unique, 0);

  BOOST_CHECK(0 == r);
}

BOOST_AUTO_TEST_CASE(Find2ByContainer)
{
  IContainer* container  = new IContainer(10,5);
  IContainer* uniqueTest = new IContainer(4,5);

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

  *uniqueTest << 1 << 2 << 3 << 4 << 5;
  *uniqueTest << 2 << 2 << 3 << 4 << 5;
  *uniqueTest << 1 << 2 << 1 << 4 << 5;
  *uniqueTest << 1 << 2 << 5 << 4 << 5;

  IContainer *unique = container->unique();
  int r = -1;

  BOOST_CHECK(uniqueTest->equals(unique));

  r = container->find(unique, 0);
  BOOST_CHECK(0 == r);

  r = container->find(unique, 1);
  BOOST_CHECK(1 == r);

  r = container->find(unique, 2);
  BOOST_CHECK(5 == r);

  r = container->find(unique, 3);
  BOOST_CHECK(9 == r);

}

BOOST_AUTO_TEST_CASE(FindList1ByContainer)
{
  IContainer* container = new IContainer(10,5);

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

  IContainer *unique = container->unique();

  vector<int> r = container->findlist(unique, 0);

  BOOST_CHECK(10 == (int)r.size());
  for(int i = 0; i < 10; i++)
  {
    BOOST_CHECK(i == r[i]);
  }
}

BOOST_AUTO_TEST_CASE(FindList2ByContainer)
{
  IContainer* container  = new IContainer(10,5);
  IContainer* uniqueTest = new IContainer(4,5);

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

  *uniqueTest << 1 << 2 << 3 << 4 << 5;
  *uniqueTest << 2 << 2 << 3 << 4 << 5;
  *uniqueTest << 1 << 2 << 1 << 4 << 5;
  *uniqueTest << 1 << 2 << 5 << 4 << 5;

  IContainer* unique = container->unique();

  BOOST_CHECK(unique->equals(uniqueTest));

  vector<int> r;

  r = container->findlist(unique, 0);
  BOOST_CHECK(5 == (int)r.size());
  BOOST_CHECK(0 == r[0]);
  BOOST_CHECK(2 == r[1]);
  BOOST_CHECK(4 == r[2]);
  BOOST_CHECK(6 == r[3]);
  BOOST_CHECK(8 == r[4]);

  r = container->findlist(unique, 1);
  BOOST_CHECK(2 == (int)r.size());
  BOOST_CHECK(1 == r[0]);
  BOOST_CHECK(3 == r[1]);

  r = container->findlist(unique, 2);
  BOOST_CHECK(2 == (int)r.size());
  BOOST_CHECK(5 == r[0]);
  BOOST_CHECK(7 == r[1]);

  r = container->findlist(unique, 3);
  BOOST_CHECK(1 == (int)r.size());
  BOOST_CHECK(9 == r[0]);
}

BOOST_AUTO_TEST_CASE(FindList3ByContainer)
{
  IContainer* container = new IContainer(10,5);

  *container << 1 << 2 << 1  << 11 << 21;
  *container << 1 << 2 << 2  << 12 << 22;
  *container << 1 << 2 << 3  << 13 << 23;
  *container << 1 << 2 << 4  << 14 << 24;
  *container << 1 << 2 << 5  << 15 << 25;
  *container << 1 << 3 << 10 << 16 << 26;
  *container << 1 << 3 << 7  << 17 << 27;
  *container << 1 << 3 << 8  << 18 << 28;
  *container << 1 << 3 << 9  << 19 << 29;
  *container << 1 << 3 << 10 << 20 << 30;

  vector<int> r;
  vector<int> indices;
  indices.push_back(0);
  indices.push_back(1);

  r = container->findlist(container, 0, indices);
  BOOST_CHECK(5 == (int)r.size());
  BOOST_CHECK(0 == r[0]);
  BOOST_CHECK(1 == r[1]);
  BOOST_CHECK(2 == r[2]);
  BOOST_CHECK(3 == r[3]);
  BOOST_CHECK(4 == r[4]);

  indices.clear();
  indices.push_back(0);
  indices.push_back(2);

  r = container->findlist(container, 5, indices);
  BOOST_CHECK(2 == (int)r.size());
  BOOST_CHECK(5 == r[0]);
  BOOST_CHECK(9 == r[1]);

}

BOOST_AUTO_TEST_CASE(FindList4ByContainer)
{
  IContainer* container = new IContainer(10,5);

  *container << 1 << 2 << 1  << 11 << 21;
  *container << 1 << 2 << 2  << 12 << 22;
  *container << 1 << 2 << 3  << 13 << 23;
  *container << 1 << 2 << 4  << 14 << 24;
  *container << 1 << 2 << 5  << 15 << 25;
  *container << 1 << 3 << 10 << 16 << 26;
  *container << 1 << 3 << 7  << 17 << 27;
  *container << 1 << 3 << 8  << 18 << 28;
  *container << 1 << 3 << 9  << 19 << 29;
  *container << 1 << 3 << 10 << 20 << 30;

  vector<int> r;
  vector<int> indices;
  indices.push_back(0);
  indices.push_back(1);
  vector<int> values;
  values.push_back(1);
  values.push_back(2);

  r = container->findlist(values, indices);
  BOOST_CHECK(5 == (int)r.size());
  BOOST_CHECK(0 == r[0]);
  BOOST_CHECK(1 == r[1]);
  BOOST_CHECK(2 == r[2]);
  BOOST_CHECK(3 == r[3]);
  BOOST_CHECK(4 == r[4]);

  indices.clear();
  indices.push_back(0);
  indices.push_back(2);
  values.clear();
  values.push_back(1);
  values.push_back(10);

  r = container->findlist(container, 5, indices);
  BOOST_CHECK(2 == (int)r.size());
  BOOST_CHECK(5 == r[0]);
  BOOST_CHECK(9 == r[1]);

}

BOOST_AUTO_TEST_CASE(GetRow1)
{
  IContainer* container = new IContainer(10,5);

  *container << 1 << 2 << 1  << 11 << 21;
  *container << 1 << 2 << 2  << 12 << 22;
  *container << 1 << 2 << 3  << 13 << 23;
  *container << 1 << 2 << 4  << 14 << 24;
  *container << 1 << 2 << 5  << 15 << 25;
  *container << 1 << 3 << 10 << 16 << 26;
  *container << 1 << 3 << 7  << 17 << 27;
  *container << 1 << 3 << 8  << 18 << 28;
  *container << 1 << 3 << 9  << 19 << 29;
  *container << 1 << 3 << 10 << 20 << 30;

  vector<int> r = container->row(0);

  BOOST_CHECK(5  == (int)r.size());
  BOOST_CHECK(1  == r[0]);
  BOOST_CHECK(2  == r[1]);
  BOOST_CHECK(1  == r[2]);
  BOOST_CHECK(11 == r[3]);
  BOOST_CHECK(21 == r[4]);

  r = container->row(1);

  BOOST_CHECK(5  == (int)r.size());
  BOOST_CHECK(1  == r[0]);
  BOOST_CHECK(2  == r[1]);
  BOOST_CHECK(2  == r[2]);
  BOOST_CHECK(12 == r[3]);
  BOOST_CHECK(22 == r[4]);

  r = container->row(2);

  BOOST_CHECK(5  == (int)r.size());
  BOOST_CHECK(1  == r[0]);
  BOOST_CHECK(2  == r[1]);
  BOOST_CHECK(3  == r[2]);
  BOOST_CHECK(13 == r[3]);
  BOOST_CHECK(23 == r[4]);

}

BOOST_AUTO_TEST_CASE(GetRow2)
{
  IContainer* container = new IContainer(10,5);

  *container << 1 << 2 << 1  << 11 << 21;
  *container << 1 << 2 << 2  << 12 << 22;
  *container << 1 << 2 << 3  << 13 << 23;
  *container << 1 << 2 << 4  << 14 << 24;
  *container << 1 << 2 << 5  << 15 << 25;
  *container << 1 << 3 << 10 << 16 << 26;
  *container << 1 << 3 << 7  << 17 << 27;
  *container << 1 << 3 << 8  << 18 << 28;
  *container << 1 << 3 << 9  << 19 << 29;
  *container << 1 << 3 << 10 << 20 << 30;

  vector<int> indices;
  indices.push_back(3);
  indices.push_back(4);
  vector<int> r = container->row(0, indices);

  BOOST_CHECK(2  == (int)r.size());
  BOOST_CHECK(11 == r[0]);
  BOOST_CHECK(21 == r[1]);

  r = container->row(1, indices);
  BOOST_CHECK(2  == (int)r.size());
  BOOST_CHECK(12 == r[0]);
  BOOST_CHECK(22 == r[1]);

  r = container->row(2, indices);
  BOOST_CHECK(2  == (int)r.size());
  BOOST_CHECK(13 == r[0]);
  BOOST_CHECK(23 == r[1]);

}

BOOST_AUTO_TEST_CASE(Equal)
{
  IContainer* A = new IContainer(10,5);
  IContainer* B = new IContainer(10,5);

  for(int r = 0; r < 10; r++)
  {
    for(int c = 0; c < 5; c++)
    {
      *A << (r + c);
      *B << (r + c);
    }
  }
  BOOST_CHECK(A->equals(B));
}

BOOST_AUTO_TEST_SUITE_END()
