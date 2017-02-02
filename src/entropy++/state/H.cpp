#include <entropy++/H.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;
using namespace entropy;

double __empericalH(ULContainer* X)
{
  assert(X->isDiscretised());

  int    maxX = 0;
  double sum  = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    int x = X->get(i,0);
    if(x > maxX) maxX = x;
  }

  maxX = maxX + 1;

  double  *px  = new double[maxX];

  for(int x = 0; x < maxX; x++) px[x] = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    int x = X->get(i, 0);
    px[x] = px[x] + 1.0;
  }

  for(int x = 0; x < maxX; x++)
  {
    px[x] = px[x] / (double)(X->rows());
  }

  sum = 0.0;
  for(int x = 0; x < maxX; x++) sum += px[x];
  assert(fabs(sum - 1.0) < 0.000001);

  double r = 0.0;
  for(int x = 0; x < maxX; x++)
  {
    if(px[x] > 0.0) r -= px[x] * log2(px[x]);
  }

  delete[] px;

  return r;
}

double H(ULContainer* X, int mode)
{
  switch(mode)
  {
    case EMPERICAL:
      return __empericalH(X);
      break;
    default:
      cerr << "H::calulate unknown mode given: " << mode << endl;
      break;
  }
  return 0.0;
}

