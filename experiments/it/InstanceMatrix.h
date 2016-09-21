#ifndef __INSTANCEMATRIX_H__
#define __INSTANCEMATRIX_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <math.h>
#include <entropy++/Matrix.h>
#include "Feature.h"
#include "ITMatrix.h"

using namespace std;

class InstanceMatrix : public ITMatrix{

public:
  //InstanceMatrix();
  InstanceMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, double valuelambda);
  ~InstanceMatrix();
  // TODO copying vectors can be expensive
  vector<int> getInstanceMatrixX(int Feati,int Featj,int delti,int deltj);
  vector<int> getInstanceMatrixY(int Feati,int Featj,int delti,int deltj);

private:
  void _getMatrix(double valuelambda);
  // TODO check for a better way to store
  vector<vector<int> > ****_mat;

};
#endif
