#ifndef __INSTANCEMATRIX_H__
#define __INSTANCEMATRIX_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <entropy++/Matrix.h>
#include <entropy++/iterativescaling/Feature.h>
#include <entropy++/iterativescaling/ITMatrix.h>

#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <math.h>

using namespace std;

class InstanceMatrix : public ITMatrix{

public:
  //InstanceMatrix();
  InstanceMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, double valuelambda);
  InstanceMatrix(ULContainer &eX, ULContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, double valuelambda);
  ~InstanceMatrix();
  // TODO copying vectors can be expensive

  vector<int> getInstanceMatrixX(int feat, int deltai, int deltaj);
  vector<int> getInstanceMatrixY(int feat, int deltai, int deltaj);

private:
  void _getMatrix(double valuelambda);
  // TODO check for a better way to store
  vector<vector<int> > ***_mat;

};
#endif
