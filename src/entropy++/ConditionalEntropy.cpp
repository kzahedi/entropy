#include <entropy++/ConditionalEntropy.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;

double __empericalH(ULContainer* X, ULContainer* Y)
{
  assert(X->isDiscretised());
  assert(Y->isDiscretised());

  int    maxX = X->max(0)+1;
  int    maxY = Y->max(0)+1;
  double sum  = 0.0;

  double  *py     = new double[maxY];
  double **pxy    = new double*[maxX];
  double **px_c_y = new double*[maxX];

  for(int x = 0; x < maxX; x++)
  {
    pxy[x]    = new double[maxY];
    px_c_y[x] = new double[maxY];
  }

  for(int y = 0; y < maxY; y++)
  {
    py[y] = 0.0;
    for(int x = 0; x < maxX; x++)
    {
      pxy[x][y]    = 0.0;
      px_c_y[x][y] = 0.0;
    }
  }

  for(int i = 0; i < Y->rows(); i++)
  {
    int x     = X->get(i, 0);
    int y     = Y->get(i, 0);
    py[y]     = py[y] + 1.0;
    pxy[x][y] = pxy[x][y] + 1.0;
  }

  for(int y = 0; y < maxY; y++)
  {
    py[y] = py[y] / (double)(Y->rows());
    for(int x = 0; x < maxX; x++)
    {
      pxy[x][y] = pxy[x][y] / (double)(X->rows());
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      px_c_y[x][y] = pxy[x][y] / py[y];
    }
  }

  sum = 0.0;
  for(int y = 0; y < maxY; y++) sum += py[y];
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int x = 0; x < maxY; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      sum += pxy[x][y];
    }
  }
  assert(fabs(sum - 1.0) < 0.000001);

  double r = 0.0;
  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      if(pxy[x][y] > 0.0 && px_c_y[x][y] > 0.0) r -= pxy[x][y] * log2(px_c_y[x][y]);
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    delete [] pxy[x];
    delete [] px_c_y[x];
  }

  delete[] py;
  delete[] pxy;
  delete[] px_c_y;

  return r;
}

double entropy::ConditionalEntropy(ULContainer* X, ULContainer* Y, int mode)
{
  switch(mode)
  {
    case EMPERICAL:
      return __empericalH(X,Y);
      break;
    default:
      cerr << "H::calulate unknown mode given: " << mode << endl;
      break;
  }
  return 0.0;
}

