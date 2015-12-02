#include <entropy++/MIs.h>

#include <entropy++/SparseMatrix.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;

MIs::MIs()
{
  _mode = EMPERICAL;
}

MIs::~MIs()
{
}

double MIs::calculate(Container* X, Container* Y)
{
  switch(_mode)
  {
    case EMPERICAL:
      return __empericalMIs(X, Y);
      break;
    default:
      cerr << "MIs::calulate unknown mode given: " << _mode << endl;
      break;
  }
  return 0.0;
}

double MIs::__empericalMIs(Container* X, Container* Y)
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

  SparseMatrix pxy;
  SparseMatrix px;
  SparseMatrix py;

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

  double r = 0.0;
  for(int i = 0; i < pxy.size(); i++)
  {
    MatrixIndex mi = pxy.getmi(i);
    int x = mi.first;
    int y = mi.second;
    if(px(x) > 0.0 && py(y) > 0.0 && pxy(x,y) > 0.0)
    {
      r += pxy(x,y) * (log2(pxy(x,y)) - log2(px(x) * py(y)));
    }
  }

  return r;
}
