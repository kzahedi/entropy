#include "entropy++/distribution/CMI.h"

#include <math.h>
#include <assert.h>

double entropy::distribution::CMI(double*** pxyz, int dimX, int dimY, int dimZ)
{
  double ***pxy_c_z = new double**[dimX];
  double  **px_c_z  = new double*[dimX];
  double  **py_c_z  = new double*[dimY];
  double   *pz      = new double[dimZ];

  for(int z = 0; z < dimZ; z++)
  {
    pz[z] = 0.0;
  }

  for(int x = 0; x < dimX; x++)
  {
    pxy_c_z[x] = new double*[dimY];
    for(int y = 0; y < dimY; y++)
    {
      pxy_c_z[x][y] = new double[dimZ];
      for(int z = 0; z < dimZ; z++)
      {
        pxy_c_z[x][y][z] = 0.0;
      }
    }
  }

  for(int x = 0; x < dimX; x++)
  {
    px_c_z[x] = new double[dimZ];
    for(int z = 0; z < dimZ; z++)
    {
      px_c_z[x][z] = 0.0;
    }
  }

  for(int y = 0; y < dimY; y++)
  {
    py_c_z[y] = new double[dimZ];
    for(int z = 0; z < dimZ; z++)
    {
      py_c_z[y][z] = 0.0;
    }
  }

  double sum = 0.0;
  for(int x = 0; x < dimX; x++)
  {
    for(int y = 0; y < dimY; y++)
    {
      for(int z = 0; z < dimZ; z++)
      {
        sum += pxyz[x][y][z];
      }
    }
  }
  assert(fabs(sum - 1.0) < 0.000001);

  for(int x = 0; x < dimX; x++)
  {
    for(int y = 0; y < dimY; y++)
    {
      for(int z = 0; z < dimZ; z++)
      {
        pz[z] = pz[z] + pxyz[x][y][z];
      }
    }
  }

  sum = 0.0;
  for(int z = 0; z < dimZ; z++) sum += pz[z];
  assert(fabs(sum - 1.0) < 0.000001);

  for(int x = 0; x < dimX; x++)
  {
    for(int y = 0; y < dimY; y++)
    {
      for(int z = 0; z < dimZ; z++)
      {
        px_c_z[x][z] = px_c_z[x][z] + pxyz[x][y][z];
        py_c_z[y][z] = py_c_z[y][z] + pxyz[x][y][z];
      }
    }
  }

  for(int x = 0; x < dimX; x++)
  {
    for(int y = 0; y < dimY; y++)
    {
      for(int z = 0; z < dimZ; z++)
      {
        if(pz[z] > 0.0)
        {
          pxy_c_z[x][y][z] = pxyz[x][y][z] / pz[z];
        }
      }
    }
  }

  for(int x = 0; x < dimX; x++)
  {
    for(int z = 0; z < dimZ; z++)
    {
      px_c_z[x][z] = px_c_z[x][z] / pz[z];
    }
  }

  for(int y = 0; y < dimY; y++)
  {
    for(int z = 0; z < dimZ; z++)
    {
      py_c_z[y][z] = py_c_z[y][z] / pz[z];
    }
  }

  for(int z = 0; z < dimZ; z++)
  {
    sum = 0.0;
    for(int x = 0; x < dimX; x++) sum += px_c_z[x][z];
    assert(fabs(sum - 1.0) < 0.000001);

    sum = 0.0;
    for(int y = 0; y < dimY; y++) sum += py_c_z[y][z];
    assert(fabs(sum - 1.0) < 0.000001);
  }

  sum = 0.0;
  for(int z = 0; z < dimZ; z++) sum += pz[z];
  assert(fabs(sum - 1.0) < 0.000001);

  double r = 0.0;
  for(int x = 0; x < dimX; x++)
  {
    for(int y = 0; y < dimY; y++)
    {
      for(int z = 0; z < dimZ; z++)
      {
        if(pxyz[x][y][z]    > 0.0 &&
           pxy_c_z[x][y][z] > 0.0 &&
           px_c_z[x][z]     > 0.0 &&
           py_c_z[y][z]     > 0.0)
        {
          r += pxyz[x][y][z]
            * (log2(pxy_c_z[x][y][z]) - log2(px_c_z[x][z] * py_c_z[y][z]));
        }
      }
    }
  }

  for(int x = 0; x < dimX; x++)
  {
    for(int y = 0; y < dimY; y++)
    {
      delete[] pxy_c_z[x][y];
    }
    delete[] pxy_c_z[x];
  }

  for(int x = 0; x < dimX; x++)
  {
    delete[] px_c_z[x];
  }

  for(int y = 0; y < dimY; y++)
  {
    delete[] py_c_z[y];
  }

  delete[] py_c_z;
  delete[] px_c_z;
  delete[] pz;

  return r;

}
