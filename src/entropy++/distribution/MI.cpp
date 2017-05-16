#include "MI.h"

#include <math.h>
#include <assert.h>

double entropy::distribution::MI(double** pxy, int dimX, int dimY)
{
  double *px = new double[dimX];
  double *py = new double[dimY];

  for(int x = 0; x < dimX; x++) px[x] = 0.0;
  for(int y = 0; y < dimY; y++) py[y] = 0.0;

  for(int x = 0; x < dimX; x++)
  {
    for(int y = 0; y < dimY; y++)
    {
      py[y] += pxy[x][y];
      px[x] += pxy[x][y];
    }
  }

  double sum = 0.0;
  for(int x = 0; x < dimX; x++) sum += px[x];
  assert(fabs(sum - 1.0) < 0.000001);

  sum = 0.0;
  for(int y = 0; y < dimY; y++) sum += py[y];
  assert(fabs(sum - 1.0) < 0.000001);

  double r = 0.0;
  for(int x = 0; x < dimX; x++)
  {
    for(int y = 0; y < dimY; y++)
    {
      if(px[x] > 0.0 && py[y] > 0.0 && pxy[x][y] > 0.0)
      {
        r += pxy[x][y] * (log2(pxy[x][y]) - log2(px[x] * py[y]));
      }
    }
  }

  delete[] px;
  delete[] py;

  return r;

}

