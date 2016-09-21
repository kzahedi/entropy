#include <entropy++/Csv.h>
#include <entropy++/sparse/MC_MI.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/sparse/state/MC_W.h>
#include <entropy++/sparse/state/MC_MI.h>
#include <entropy++/Matrix.h>
#include <entropy++/SparseMatrix.h>

#include <time.h>
#include <string>
#include <iostream>
#include <math.h>
#include <vector>
#include "Feature.h"
#include "FeatureMatrix.h"
#include "GIS.h"
#include "InstanceMatrix.h"
#include "FeatureMatrixsp.h"
#include "SCGIS.h"
#include "GISsp.h"
#include "Test.h"



int main(int argc, char **argv)
{
  DContainer *zX = new DContainer(2,1); // alphabet
  *zX << 0 << 1;
  DContainer *zY = new DContainer(2,1); // alphabet
  *zY << 0 << 1;
  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;

  ItParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;

  param.maxit          = 500;
  param.konv           = 0.001;
  param.time           = false;
  param.test           = true;
  param.seconds        = 1;

  Test *test = new Test(2, 2, 100, lambda, *zX, *zY, param); // for test cases
  // Test *test = new Test(5, 5, 100000, lambda, *zX, *zY, 500, 0.0001, true, 2);
  test->comparison();
}
