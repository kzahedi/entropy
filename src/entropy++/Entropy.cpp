#include <entropy++/Entropy.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;

Entropy::Entropy()
{
  _mode = Entropy_EMPERICAL;
}

Entropy::~Entropy()
{
}

double Entropy::calulate(Container* X)
{
  switch(_mode)
  {
    case Entropy_EMPERICAL:
      return __empericalEntropy(X);
      break;
    default:
      cerr << "Entropy::calulate unknown mode given: " << _mode << endl;
      break;
  }
  return 0.0;
}

double Entropy::__empericalEntropy(Container* X)
{
  assert(X->isDiscretised());

  int    maxX = 0;
  double sum  = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    if(X->get(i,0) > maxX) maxX = X->get(i,0);
  }

  maxX = maxX + 1;

  double *px = new double[maxX];

  for(int x = 0; x < maxX; x++) px[x] = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    int x = X->get(i, 0);
    px[x] = px[x] + 1.0;
  }
  sum = 0.0;
  for(int x = 0; x < maxX; x++)
  {
    sum += px[x];
  }
  assert(fabs(sum - X->rows()) < 0.000001);

  for(int x = 0; x < maxX; x++)
  {
    px[x] = px[x] / (double)(X->rows());
  }


  double r = 0.0;
  for(int x = 0; x < maxX; x++)
  {
    if(px[x] > 0.0) r -= px[x] * log2(px[x]);
  }

  delete[] px;

  return r;
}
