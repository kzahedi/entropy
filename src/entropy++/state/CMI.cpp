#include <entropy++/state/CMI.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;
using namespace entropy;
using namespace entropy::state;


DContainer* __empericalCMIsd(ULContainer* X, ULContainer* Y, ULContainer* Z)
{
  assert(X->isDiscretised());
  assert(Y->isDiscretised());
  assert(Z->isDiscretised());

  int    maxX = 0;
  int    maxY = 0;
  int    maxZ = 0;
  double sum  = 0.0;

  for(int i = 0; i < X->rows(); i++)
  {
    if(X->get(i,0) > maxX) maxX = X->get(i,0);
    if(Y->get(i,0) > maxY) maxY = Y->get(i,0);
    if(Z->get(i,0) > maxZ) maxZ = Z->get(i,0);
  }

  maxX = maxX + 1;
  maxY = maxY + 1;
  maxZ = maxZ + 1;

  double ***pxyz    = new double**[maxX];
  double ***pxy_c_z = new double**[maxX];
  double  **px_c_z  = new double*[maxX];
  double  **py_c_z  = new double*[maxY];
  double   *pz      = new double[maxZ];
  double ***cmi     = new double**[maxX];

  for(int z = 0; z < maxZ; z++)
  {
    pz[z] = 0.0;
  }

  for(int x = 0; x < maxX; x++)
  {
    pxyz[x]    = new double*[maxY];
    cmi[x]     = new double*[maxY];
    pxy_c_z[x] = new double*[maxY];
    for(int y = 0; y < maxY; y++)
    {
      pxyz[x][y]    = new double[maxZ];
      cmi[x][y]     = new double[maxZ];
      pxy_c_z[x][y] = new double[maxZ];
      for(int z = 0; z < maxZ; z++)
      {
        pxyz[x][y][z]    = 0.0;
        cmi[x][y][z]     = 0.0;
        pxy_c_z[x][y][z] = 0.0;
      }
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    px_c_z[x] = new double[maxZ];
    for(int z = 0; z < maxZ; z++)
    {
      px_c_z[x][z] = 0.0;
    }
  }

  for(int y = 0; y < maxY; y++)
  {
    py_c_z[y] = new double[maxZ];
    for(int z = 0; z < maxZ; z++)
    {
      py_c_z[y][z] = 0.0;
    }
  }

  for(int i = 0; i < X->rows(); i++)
  {
    int x         = X->get(i, 0);
    int y         = Y->get(i, 0);
    int z         = Z->get(i, 0);
    pxyz[x][y][z] = pxyz[x][y][z] + 1.0;
  }

  sum = 0.0;
  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      for(int z = 0; z < maxZ; z++)
      {
        pxyz[x][y][z] = pxyz[x][y][z] / (double)(X->rows());
        sum += pxyz[x][y][z];
      }
    }
  }
  assert(fabs(sum - 1.0) < 0.000001);

  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      for(int z = 0; z < maxZ; z++)
      {
        pz[z] = pz[z] + pxyz[x][y][z];
      }
    }
  }

  sum = 0.0;
  for(int z = 0; z < maxZ; z++) sum += pz[z];
  assert(fabs(sum - 1.0) < 0.000001);

  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      for(int z = 0; z < maxZ; z++)
      {
        px_c_z[x][z] = px_c_z[x][z] + pxyz[x][y][z];
        py_c_z[y][z] = py_c_z[y][z] + pxyz[x][y][z];
      }
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      for(int z = 0; z < maxZ; z++)
      {
        if(pz[z] > 0.0)
        {
          pxy_c_z[x][y][z] = pxyz[x][y][z] / pz[z];
        }
      }
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    for(int z = 0; z < maxZ; z++)
    {
      px_c_z[x][z] = px_c_z[x][z] / pz[z];
    }
  }

  for(int y = 0; y < maxY; y++)
  {
    for(int z = 0; z < maxZ; z++)
    {
      py_c_z[y][z] = py_c_z[y][z] / pz[z];
    }
  }

  for(int z = 0; z < maxZ; z++)
  {
    sum = 0.0;
    for(int x = 0; x < maxX; x++) sum += px_c_z[x][z];
    assert(fabs(sum - 1.0) < 0.000001);

    sum = 0.0;
    for(int y = 0; y < maxY; y++) sum += py_c_z[y][z];
    assert(fabs(sum - 1.0) < 0.000001);
  }

  sum = 0.0;
  for(int z = 0; z < maxZ; z++) sum += pz[z];
  assert(fabs(sum - 1.0) < 0.000001);

  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      for(int z = 0; z < maxZ; z++)
      {
        if(pxyz[x][y][z]    > 0.0 &&
           pxy_c_z[x][y][z] > 0.0 &&
           px_c_z[x][z]     > 0.0 &&
           py_c_z[y][z]     > 0.0)
        {
          cmi[x][y][z] = (log2(pxy_c_z[x][y][z]) - log2(px_c_z[x][z] * py_c_z[y][z]));
        }
      }
    }
  }

  DContainer* r = new DContainer(X->rows(), 1);
  for(int i = 0; i < X->rows(); i++)
  {
    int x     = X->get(i, 0);
    int y     = Y->get(i, 0);
    int z     = Z->get(i, 0);
    (*r)(i,0) = cmi[x][y][z];
  }

  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      delete[] pxyz[x][y];
      delete[] cmi[x][y];
      delete[] pxy_c_z[x][y];
    }
    delete[] pxyz[x];
    delete[] cmi[x];
    delete[] pxy_c_z[x];
  }

  for(int x = 0; x < maxX; x++)
  {
    delete[] px_c_z[x];
  }

  for(int y = 0; y < maxY; y++)
  {
    delete[] py_c_z[y];
  }

  delete[] pxyz;
  delete[] cmi;
  delete[] py_c_z;
  delete[] px_c_z;
  delete[] pz;

  return r;
}

DContainer* entropy::state::CMI(ULContainer* X, ULContainer* Y, ULContainer *Z, int mode)
{
  switch(mode)
  {
    case EMPERICAL:
      return __empericalCMIsd(X, Y, Z);
      break;
    default:
      cerr << "CMIsd::calulate unknown mode given: " << mode << endl;
      break;
  }
  return NULL;
}

