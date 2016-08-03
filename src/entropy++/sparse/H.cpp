#include <entropy++/sparse/H.h>

#include <entropy++/SparseMatrix.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;
using namespace entropy::sparse;

double __empericalHs(ULContainer* X)
{
  assert(X->isDiscretised());

  int    maxX = 0;
  double sum  = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    if(X->get(i,0) > maxX) maxX = X->get(i,0);
  }

  maxX = maxX + 1;

  SparseMatrix px;

  for(int i = 0; i < X->rows(); i++)
  {
    int x    = X->get(i, 0);
    px(x)    = px(x)    + 1.0;
  }

  px /= (double)(X->rows());

  sum = 0.0;
  for(int i = 0; i < (int)px.size(); i++) sum += px.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  double r = 0.0;
  for(int i = 0; i < px.size(); i++)
  {
    MatrixIndex mi = px.getmi(i);
    int x = mi.first;
    // cout << "px(" << x << ") = " << px(x) << endl;
    if(px(x) > 0.0) r += px(x) * log2(px(x));
  }

  return -r;
}

double entropy::sparse::H(ULContainer* X, int mode)
{
  switch(mode)
  {
    case EMPERICAL:
      return __empericalHs(X);
      break;
    default:
      cerr << "Hs::calulate unknown mode given: " << mode << endl;
      break;
  }
  return 0.0;
}

