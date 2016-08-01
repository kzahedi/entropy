#include <entropy++/sparse/CMI.h>
#include <entropy++/SparseMatrix.h>

#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;
using namespace entropy::sparse;


double __empericalCMIs(Container* X, Container* Y, Container* Z)
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

  SparseMatrix pxyz;
  SparseMatrix pxy_c_z;
  SparseMatrix px_c_z;
  SparseMatrix py_c_z;
  SparseMatrix pz;

  for(int i = 0; i < X->rows(); i++)
  {
    int x = X->get(i, 0);
    int y = Y->get(i, 0);
    int z = Z->get(i, 0);
    pxyz(x,y,z) = pxyz(x,y,z) + 1.0;
  }

  pxyz /= (double)(X->rows());

  sum = 0.0;
  for(int i = 0; i < pxyz.size(); i++) sum += pxyz.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  for(int i = 0; i < pxyz.size(); i++)
  {
    MatrixIndex mi = pxyz.getmi(i);
    int x = mi.first;
    int y = mi.second;
    int z = mi.third;
    pz(z) += pxyz(x,y,z);
  }

  sum = 0.0;
  for(int i = 0; i < pz.size(); i++) sum += pz.get(i);
  assert(fabs(sum - 1.0) < 0.000001);

  for(int i = 0; i < pxyz.size(); i++)
  {
    MatrixIndex mi = pxyz.getmi(i);
    int x = mi.first;
    int y = mi.second;
    int z = mi.third;
    if(pz(z) > 0.0) pxy_c_z(x,y,z) = pxyz(x,y,z) / pz(z);
  }

  for(int i = 0; i < pxyz.size(); i++)
  {
    MatrixIndex mi = pxyz.getmi(i);
    int x = mi.first;
    int y = mi.second;
    int z = mi.third;
    py_c_z(y,z) += pxyz(x,y,z);
  }
  for(int i = 0; i < py_c_z.size(); i++)
  {
    MatrixIndex mi = py_c_z.getmi(i);
    int y = mi.first;
    int z = mi.second;
    py_c_z(y,z) = py_c_z(y,z) / pz(z);
  }

  for(int i = 0; i < pxyz.size(); i++)
  {
    MatrixIndex mi = pxyz.getmi(i);
    int x = mi.first;
    int y = mi.second;
    int z = mi.third;
    px_c_z(x,z) += pxyz(x,y,z);
  }
  for(int i = 0; i < px_c_z.size(); i++)
  {
    MatrixIndex mi = px_c_z.getmi(i);
    int x = mi.first;
    int z = mi.second;
    px_c_z(x,z) = px_c_z(x,z) / pz(z);
  }


  // CHECK MATRICES

  double r = 0.0;
  for(int i = 0; i < pxyz.size(); i++)
  {
    MatrixIndex mi = pxyz.getmi(i);
    int x = mi.first;
    int y = mi.second;
    int z = mi.third;
    if(pxyz(x,y,z) > 0.0 && pxy_c_z(x,y,z) > 0.0 &&
       px_c_z(x,z) > 0.0 && py_c_z(y,z)    > 0.0)
    {
      r += pxyz(x,y,z) * (log2(pxy_c_z(x,y,z)) - log2(px_c_z(x,z) * py_c_z(y,z)));
    }
  }

  return r;
}

//
// I(X;Y|Z) = \sum_{x,y,z} p(x,y,z) log( p(x,y|z) / (p(x|z) * p(y|z)))
//
double entropy::sparse::CMI(Container* X, Container* Y, Container *Z, int mode)
{
  switch(mode)
  {
    case EMPERICAL:
      return __empericalCMIs(X, Y, Z);
      break;
    default:
      cerr << "CMIs::calulate unknown mode given: " << mode << endl;
      break;
  }
  return 0.0;
}