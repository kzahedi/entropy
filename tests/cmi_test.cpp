#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cmi_test
#include <boost/test/unit_test.hpp>

#include <entropy++/Container.h>
#include <entropy++/CMI.h>
#include <entropy++/sparse/CMI.h>
#include <entropy++/SparseMatrix.h>

#include <iostream>
#include <string>

#include <math.h>

using namespace std;
using namespace entropy;


BOOST_AUTO_TEST_CASE(Sinus)
{

  DContainer X(1000,1);
  DContainer Y(1000,1);
  DContainer Z(1000,1);

  for(float i = 0; i < 1000.0; i = i + 1.0)
  {
    X << cos(i/10.0);
    Y << sin(i/5.0);
    Z << cos(i/5.0) * sin(i/5.0);
  }

  double **dom = new double*[1];
  dom[0]       = new double[2];
  dom[0][0]    = -1.0;
  dom[0][1]    =  1.0;

  int *bins    = new int[1];
  bins[0]      = 100;

  X.setDomains(dom);
  X.setBinSizes(bins);

  Y.setDomains(dom);
  Y.setBinSizes(bins);

  Z.setDomains(dom);
  Z.setBinSizes(bins);

  ULContainer *dx = X.discretise();
  ULContainer *dy = Y.discretise();
  ULContainer *dz = Z.discretise();

  double s = entropy::CMI(dx, dy, dz);

  BOOST_CHECK_CLOSE(2.473493, s, 0.001); // recalcuate somewhere else

  delete dx;
  delete dy;
}


BOOST_AUTO_TEST_CASE(SparseVsNonSparse)
{
  DContainer X(1000,1);
  DContainer Y(1000,1);
  DContainer Z(1000,1);

  for(float i = 0; i < 1000.0; i = i + 1.0)
  {
    X << cos(i/10.0);
    Y << sin(i/5.0);
    Z << cos(i/5.0) * sin(i/5.0);
  }

  double **dom = new double*[1];
  dom[0]       = new double[2];
  dom[0][0]    = -1.0;
  dom[0][1]    =  1.0;

  int *bins    = new int[1];
  bins[0]      = 100;

  X.setDomains(dom);
  X.setBinSizes(bins);

  Y.setDomains(dom);
  Y.setBinSizes(bins);

  Z.setDomains(dom);
  Z.setBinSizes(bins);

  ULContainer *dx = X.discretise();
  ULContainer *dy = Y.discretise();
  ULContainer *dz = Z.discretise();

  double s1 = entropy::CMI(dx, dy, dz);
  double s2 = entropy::sparse::CMI(dx, dy, dz);
  BOOST_CHECK_CLOSE(s1, s2, 0.001);

  delete dx;
  delete dy;
}


BOOST_AUTO_TEST_CASE(MatrixWiseComparision)
{
  //
  // Generating data
  //
  DContainer X(1000,1);
  DContainer Y(1000,1);
  DContainer Z(1000,1);

  for(float i = 0; i < 1000.0; i = i + 1.0)
  {
    X << cos(i/10.0);
    Y << sin(i/5.0);
    Z << cos(i/5.0) * sin(i/5.0);
  }

  double **dom = new double*[1];
  dom[0]       = new double[2];
  dom[0][0]    = -1.0;
  dom[0][1]    =  1.0;

  int *bins    = new int[1];
  bins[0]      = 10;

  X.setDomains(dom);
  X.setBinSizes(bins);

  Y.setDomains(dom);
  Y.setBinSizes(bins);

  Z.setDomains(dom);
  Z.setBinSizes(bins);

  ULContainer *dx = X.discretise();
  ULContainer *dy = Y.discretise();
  ULContainer *dz = Z.discretise();

  //
  // array version
  //
  int maxX = 0;
  int maxY = 0;
  int maxZ = 0;

  for(int i = 0; i < X.rows(); i++)
  {
    if(dx->get(i,0) > maxX) maxX = dx->get(i,0);
    if(dy->get(i,0) > maxY) maxY = dy->get(i,0);
    if(dz->get(i,0) > maxZ) maxZ = dz->get(i,0);
  }

  maxX = maxX + 1;
  maxY = maxY + 1;
  maxZ = maxZ + 1;

  double ***pxyz    = new double**[maxX];
  double ***pxy_c_z = new double**[maxX];
  double  **px_c_z  = new double*[maxX];
  double  **py_c_z  = new double*[maxY];
  double   *pz      = new double[maxZ];

  for(int z = 0; z < maxZ; z++)
  {
    pz[z] = 0.0;
  }

  for(int x = 0; x < maxX; x++)
  {
    pxyz[x]    = new double*[maxY];
    pxy_c_z[x] = new double*[maxY];
    for(int y = 0; y < maxY; y++)
    {
      pxyz[x][y]    = new double[maxZ];
      pxy_c_z[x][y] = new double[maxZ];
      for(int z = 0; z < maxZ; z++)
      {
        pxyz[x][y][z]    = 0.0;
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


  // Test matrices
  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      for(int z = 0; z < maxZ; z++)
      {
        BOOST_CHECK_CLOSE(0.0, pxyz[x][y][z], 0.001);
        BOOST_CHECK_CLOSE(0.0, pxy_c_z[x][y][z], 0.001);
      }
    }
  }

  for(int y = 0; y < maxY; y++)
  {
    for(int z = 0; z < maxZ; z++)
    {
      BOOST_CHECK_CLOSE(0.0, py_c_z[y][z], 0.001);
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    for(int z = 0; z < maxZ; z++)
    {
      BOOST_CHECK_CLOSE(0.0, px_c_z[x][z], 0.001);
    }
  }
  for(int z = 0; z < maxZ; z++)
  {
    BOOST_CHECK_CLOSE(0.0, pz[z], 0.001);
  }

  for(int i = 0; i < dx->rows(); i++)
  {
    int x         = dx->get(i, 0);
    int y         = dy->get(i, 0);
    int z         = dz->get(i, 0);
    pxyz[x][y][z] = pxyz[x][y][z] + 1.0;
  }

  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      for(int z = 0; z < maxZ; z++)
      {
        pxyz[x][y][z] = pxyz[x][y][z] / (double)(X.rows());
      }
    }
  }

  double sum = 0.0;
  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      for(int z = 0; z < maxZ; z++)
      {
        sum += pxyz[x][y][z];
        BOOST_TEST((pxyz[x][y][z] >= 0.0 && pxyz[x][y][z] <= 1.0));
      }
    }
  }
  BOOST_CHECK_CLOSE(1.0, sum, 0.001);

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
        if(pz[z] > 0.0) pxy_c_z[x][y][z] = pxyz[x][y][z] / pz[z];
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

  double array_r = 0.0;
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
          array_r += pxyz[x][y][z]
            * (log2(pxy_c_z[x][y][z]) - log2(px_c_z[x][z] * py_c_z[y][z]));
        }
      }
    }
  }


  //
  // Sparse Matrix version
  //

  entropy::SparseMatrix spxyz;
  entropy::SparseMatrix spxy_c_z;
  entropy::SparseMatrix spx_c_z;
  entropy::SparseMatrix spy_c_z;
  entropy::SparseMatrix spz;

  for(int i = 0; i < dx->rows(); i++)
  {
    int x = dx->get(i, 0);
    int y = dy->get(i, 0);
    int z = dz->get(i, 0);
    spxyz(x,y,z) = spxyz(x,y,z) + 1.0;
    spz(z)       = spz(z)       + 1.0;
  }

  spxyz /= (double)(X.rows());
  spz   /= (double)(X.rows());

  for(int i = 0; i < spxyz.size(); i++)
  {
    entropy::MatrixIndex mi = spxyz.getmi(i);
    int x = mi.first;
    int y = mi.second;
    int z = mi.third;
    if(spz(z) > 0.0) spxy_c_z(x,y,z) = spxyz(x,y,z) / spz(z);
  }

  for(int i = 0; i < spxyz.size(); i++)
  {
    entropy::MatrixIndex mi = spxyz.getmi(i);
    int x = mi.first;
    int y = mi.second;
    int z = mi.third;
    spy_c_z(y,z) = spy_c_z(y,z) + spxyz(x,y,z);
  }
  for(int i = 0; i < spy_c_z.size(); i++)
  {
    entropy::MatrixIndex mi = spy_c_z.getmi(i);
    int y = mi.first;
    int z = mi.second;
    spy_c_z(y,z) = spy_c_z(y,z) / spz(z);
  }

  for(int i = 0; i < spxyz.size(); i++)
  {
    entropy::MatrixIndex mi = spxyz.getmi(i);
    int x = mi.first;
    int y = mi.second;
    int z = mi.third;
    spx_c_z(x,z) = spx_c_z(x,z) + spxyz(x,y,z);
  }
  for(int i = 0; i < spx_c_z.size(); i++)
  {
    entropy::MatrixIndex mi = spx_c_z.getmi(i);
    int x = mi.first;
    int z = mi.second;
    spx_c_z(x,z) = spx_c_z(x,z) / spz(z);
  }


  // CHECK MATRICES

  double s_r = 0.0;
  for(int i = 0; i < spxyz.size(); i++)
  {
    entropy::MatrixIndex mi = spxyz.getmi(i);
    int x = mi.first;
    int y = mi.second;
    int z = mi.third;
    if(spxyz(x,y,z) > 0.0 && spxy_c_z(x,y,z) > 0.0 &&
       spx_c_z(x,z) > 0.0 && spy_c_z(y,z)    > 0.0)
    {
      s_r += spxyz(x,y,z)
             * (log2(spxy_c_z(x,y,z)) - log2(spx_c_z(x,z) * spy_c_z(y,z)));
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    for(int y = 0; y < maxY; y++)
    {
      for(int z = 0; z < maxZ; z++)
      {
        BOOST_CHECK_CLOSE(pxyz[x][y][z],    spxyz(x,y,z),    0.001);
        BOOST_CHECK_CLOSE(pxyz[x][y][z],    spxyz(x,y,z),    0.001);
        BOOST_CHECK_CLOSE(pxy_c_z[x][y][z], spxy_c_z(x,y,z), 0.001);
      }
    }
  }

  for(int x = 0; x < maxX; x++)
  {
    for(int z = 0; z < maxZ; z++)
    {
      BOOST_CHECK_CLOSE(px_c_z[x][z], spx_c_z(x,z), 0.001);
    }
  }
  for(int y = 0; y < maxY; y++)
  {
    for(int z = 0; z < maxZ; z++)
    {
      BOOST_CHECK_CLOSE(py_c_z[y][z], spy_c_z(y,z), 0.001);
    }
  }

  BOOST_CHECK_CLOSE(array_r, s_r, 0.001);

}
