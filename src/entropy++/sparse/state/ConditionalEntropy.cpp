#include <entropy++/sparse/state/ConditionalEntropy.h>

#include <entropy++/SparseMatrix.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;
using namespace entropy;

DContainer* __empericalHss(ULContainer* X, ULContainer* Y)
{
  assert(X->isDiscretised());
  assert(Y->isDiscretised());

  int    maxX = X->max(0) + 1;
  int    maxY = Y->max(0) + 1;
  double sum  = 0.0;

  entropy::SparseMatrix  pxy;
  entropy::SparseMatrix  py;
  for(int i = 0; i < X->rows(); i++)
  {
    int x = X->get(i, 0);
    int y = Y->get(i, 0);
    pxy(x,y) = pxy(x,y) + 1.0;
    py(y) = py(y) + 1.0;
  }

  pxy /= (double)(X->rows());
  py  /= (double)(Y->rows());

  sum = 0.0;
  for(int i = 0; i < (int)pxy.size(); i++) sum += pxy.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int i = 0; i < (int)py.size(); i++) sum += py.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  entropy::SparseMatrix  px_c_y;
  for(int i = 0; i < (int)pxy.size(); i++)
  {
    entropy::MatrixIndex mi = pxy.getmi(i);
    int x = mi.first;
    int y = mi.second;
    px_c_y(x,y) = pxy(x,y) / py(y);
  }

  DContainer *r = new DContainer(X->rows(), 1);

  for(int i = 0; i < X->rows(); i++)
  {
    int x = X->get(i,0);
    int y = Y->get(i,0);
    if(pxy(x,y) > 0.0 && px_c_y(x,y)) (*r) << - log2(px_c_y(x,y));
  }

  return r;
}

DContainer* entropy::sparse::state::ConditionalEntropy(ULContainer* X, ULContainer* Y, int mode)
{
  switch(mode)
  {
    case EMPERICAL:
      return __empericalHss(X,Y);
      break;
    default:
      cerr << "ConditionalEntropy unknown mode given: " << mode << endl;
      break;
  }
  return NULL;
}

