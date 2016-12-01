#include <entropy++/sparse/state/MI.h>

#include <entropy++/SparseMatrix.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;
using namespace entropy;
using namespace entropy::sparse::state;


DContainer* __empericalMIssd(ULContainer* X, ULContainer* Y)
{
  assert(X->isDiscretised());
  assert(Y->isDiscretised());

  int    maxX = 0;
  int    maxY = 0;
  double sum  = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    if(X->get(i,0) > maxX) maxX = X->get(i,0);
    if(Y->get(i,0) > maxY) maxY = Y->get(i,0);
  }

  maxX = maxX + 1;
  maxY = maxY + 1;

  entropy::SparseMatrix pxy;
  entropy::SparseMatrix mi;
  entropy::SparseMatrix px;
  entropy::SparseMatrix py;

  for(int i = 0; i < X->rows(); i++)
  {
    int x    = X->get(i, 0);
    int y    = Y->get(i, 0);
    pxy(x,y) = pxy(x,y) + 1.0;
    px(x)    = px(x)    + 1.0;
    py(y)    = py(y)    + 1.0;
  }

  pxy /= (double)(X->rows());
  px  /= (double)(X->rows());
  py  /= (double)(X->rows());

  sum = 0.0;
  for(int i = 0; i < (int)pxy.size(); i++) sum += pxy.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int i = 0; i < (int)px.size(); i++) sum += px.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int i = 0; i < (int)py.size(); i++) sum += py.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  for(int i = 0; i < pxy.size(); i++)
  {
    entropy::MatrixIndex m = pxy.getmi(i);
    int x = m.first;
    int y = m.second;
    if(px(x) > 0.0 && py(y) > 0.0 && pxy(x,y) > 0.0)
    {
      mi(x,y) = (log2(pxy(x,y)) - log2(px(x) * py(y)));
    }
  }

  DContainer *r = new DContainer(X->rows(), 1);

  for(int i = 0; i < X->rows(); i++)
  {
    int x    = X->get(i, 0);
    int y    = Y->get(i, 0);
    double v = mi(x,y);
    (*r)(i,0) = v;
  }

  return r;
}

DContainer* entropy::sparse::state::MI(ULContainer* X, ULContainer* Y, int mode)
{
  switch(mode)
  {
    case EMPERICAL:
      return __empericalMIssd(X, Y);
      break;
    default:
      cerr << "MIssd::calulate unknown mode given: " << mode << endl;
      break;
  }
  return NULL;
}

