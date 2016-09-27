#include <entropy++/PI.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;

// PI::PI()
// {
  // _mode = EMPERICAL;
// }

// PI::~PI()
// {
// }


double __empericalPI(ULContainer* X)
{
  assert(X->isDiscretised());

  int    maxX = 0;
  double sum  = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    if(X->get(i,0) > maxX) maxX = X->get(i,0);
  }

  maxX = maxX + 1;

  double **pxxp = new double*[maxX];
  double  *px   = new double[maxX];
  double  *pxp  = new double[maxX];

  for(int x = 0; x < maxX; x++)
  {
    pxxp[x] = new double[maxX];
    for(int xp = 0; xp < maxX; xp++)
    {
      pxxp[x][xp] = 0.0;
    }
  }

  for(int i = 0; i < X->rows()-1; i++)
  {
    int x       = X->get(i,     0);
    int xp      = X->get(i + 1, 0);
    pxxp[x][xp] = pxxp[x][xp] + 1.0;
    px[x]       = px[x]       + 1.0;
    pxp[xp]     = pxp[xp]     + 1.0;
  }

  for(int x = 0; x < maxX; x++)
  {
    for(int xp = 0; xp < maxX; xp++)
    {
      pxxp[x][xp] = pxxp[x][xp] / (double)(X->rows()-1);
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    px[x] = px[x] / (double)(X->rows()-1);
  }

  for(int xp = 0; xp < maxX; xp++)
  {
    pxp[xp] = pxp[xp] / (double)(X->rows()-1);
  }


  sum = 0.0;
  for(int x = 0; x < maxX; x++)
  {
    for(int xp = 0; xp < maxX; xp++)
    {
      sum += pxxp[x][xp];
    }
  }
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int x = 0; x < maxX; x++) sum += px[x];
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int xp = 0; xp < maxX; xp++) sum += pxp[xp];
  assert(fabs(sum - 1.0) < 0.000001);

  double r = 0.0;
  for(int x = 0; x < maxX; x++)
  {
    for(int xp = 0; xp < maxX; xp++)
    {
      if(px[x] > 0.0 && pxp[xp] > 0.0 && pxxp[x][xp] > 0.0)
      {
        r += pxxp[x][xp] * (log2(pxxp[x][xp]) - log2(px[x] * pxp[xp]));
      }
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    delete[] pxxp[x];
  }
  delete[] pxxp;

  return r;
}

double entropy::PI(ULContainer* X, int mode)
{
  switch(mode)
  {
    case EMPERICAL:
      return __empericalPI(X);
      break;
    default:
      cerr << "PI::calulate unknown mode given: " << mode << endl;
      break;
  }
  return 0.0;
}

