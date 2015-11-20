#include <entropy++/PIs.h>
#include <entropy++/SparseMatrix.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;

PIs::PIs()
{
  _mode = PIs_EMPERICAL;
}

PIs::~PIs()
{
}

double PIs::calculate(Container* X)
{
  switch(_mode)
  {
    case PIs_EMPERICAL:
      return __empericalPIs(X);
      break;
    default:
      cerr << "PIs::calulate unknown mode given: " << _mode << endl;
      break;
  }
  return 0.0;
}

double PIs::__empericalPIs(Container* X)
{
  assert(X->isDiscretised());

  int    maxX = 0;
  double sum  = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    if(X->get(i,0) > maxX) maxX = X->get(i,0);
  }

  maxX = maxX + 1;

  SparseMatrix pxxp;
  SparseMatrix px;
  SparseMatrix pxp;

  // double **pxxp = new double*[maxX];
  // double  *px   = new double[maxX];
  // double  *pxp  = new double[maxX];

  for(int i = 0; i < X->rows()-1; i++)
  {
    int x      = X->get(i,     0);
    int xp     = X->get(i + 1, 0);
    pxxp(x,xp) = pxxp(x,xp) + 1.0;
    px(x)      = px(x)      + 1.0;
    pxp(xp)    = pxp(xp)    + 1.0;
  }

  pxxp /= (double)(X->rows() - 1);
  px   /= (double)(X->rows() - 1);
  pxp  /= (double)(X->rows() - 1);

  sum = 0.0;
  for(int i = 0; i < (int)pxxp.size(); i++) sum += pxxp.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int i = 0; i < (int)px.size(); i++) sum += px.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int i = 0; i < (int)pxp.size(); i++) sum += pxp.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  double r = 0.0;
  for(int i = 0; i < pxxp.size(); i++)
  {
    MatrixIndex mi = pxxp.getmi(i);
    int x  = mi.first;
    int xp = mi.second;
    if(px(x) > 0.0 && pxp(xp) > 0.0 && pxxp(x,xp) > 0.0)
    {
      r += pxxp(x,xp) * (log2(pxxp(x,xp)) - log2(px(x) * pxp(xp)));
    }
  }

  return r;
}
