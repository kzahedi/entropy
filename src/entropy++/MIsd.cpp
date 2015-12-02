#include <entropy++/MIsd.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;

MIsd::MIsd()
{
  _mode = MI_EMPERICAL;
}

MIsd::~MIsd()
{
}

Container* MIsd::calculate(Container* X, Container* Y)
{
  switch(_mode)
  {
    case MI_EMPERICAL:
      return __empericalMI(X, Y);
      break;
    default:
      cerr << "MIsd::calulate unknown mode given: " << _mode << endl;
      break;
  }
  return NULL;
}

Container* MIsd::__empericalMI(Container* X, Container* Y)
{
  assert(X->isDiscretised());
  assert(Y->isDiscretised());

  int    maxX = 0;
  int    maxY = 0;
  double sum  = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    int x = X->get(i,0);
    int y = Y->get(i,0);
    if(x > maxX) maxX = x;
    if(y > maxY) maxY = y;
  }

  maxX = maxX + 1;
  maxY = maxY + 1;

  double **pxy  = new double*[maxX];
  double  *px   = new double[maxX];
  double  *py   = new double[maxY];
  double **misd = new double*[maxX];

  for(int x = 0; x < maxX; x++)
  {
    pxy[x] = new double[maxY];
    for(int y = 0; y < maxY; y++)
    {
      pxy[x][y] = 0.0;
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    misd[x] = new double[maxY];
    for(int y = 0; y < maxY; y++)
    {
      misd[x][y] = 0.0;
    }
  }

  for(int x = 0; x < maxX; x++) px[x] = 0.0;
  for(int y = 0; y < maxY; y++) py[y] = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    int x     = X->get(i, 0);
    int y     = Y->get(i, 0);
    pxy[x][y] = pxy[x][y] + 1.0;
    px[x]     = px[x]     + 1.0;
    py[y]     = py[y]     + 1.0;
  }

  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      pxy[x][y] = pxy[x][y] / (double)(X->rows());
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    px[x] = px[x] / (double)(X->rows());
  }

  for(int y = 0; y < maxY; y++)
  {
    py[y] = py[y] / (double)(Y->rows());
  }

  sum = 0.0;
  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      sum += pxy[x][y];
    }
  }
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int x = 0; x < maxX; x++) sum += px[x];
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int y = 0; y < maxY; y++) sum += py[y];
  assert(fabs(sum - 1.0) < 0.000001);

  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      if(px[x] > 0.0 && py[y] > 0.0 && pxy[x][y] > 0.0)
      {
        misd[x][y] = log2(pxy[x][y]) - log2(px[x] * py[y]);
      }
    }
  }

  Container* r_misd = new Container(X->rows(), 1);
  for(int i = 0; i < X->rows(); i++)
  {
    int x = X->get(i,0);
    int y = Y->get(i,0);
    (*r_misd)(i,0) = misd[x][y];
  }

  for(int x = 0; x < maxX; x++)
  {
    delete[] pxy[x];
    delete[] misd[x];
  }
  delete[] pxy;
  delete[] misd;

  return r_misd;
}
