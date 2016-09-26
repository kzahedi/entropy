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
  InstanceMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, double valuelambda);
  ~InstanceMatrix();
  // TODO copying vectors can be expensive

  vector<int> getInstanceMatrixDeltaX(int feat);
  vector<int> getInstanceMatrixDeltaY(int feat);
  vector<int> getInstanceMatrixX(int feat);
  vector<int> getInstanceMatrixY(int feat);

private:
  void _getMatrix(double valuelambda);
  // TODO check for a better way to store
  vector<vector<int> > *_mat;

};
#endif
