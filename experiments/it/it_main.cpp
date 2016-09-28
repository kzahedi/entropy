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
	//V(3,vector<int>(0))
vector<vector<int > > alphX(3,vector<int>(0));
alphX[0].push_back(0);
alphX[1].push_back(1);
alphX[2].push_back(2);
alphX[3].push_back(3);

vector<vector<int > > alphY(3,vector<int>(0));
alphY[0].push_back(0);
alphY[1].push_back(1);
alphY[2].push_back(2);
alphY[3].push_back(3);
  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 10;
  param.konv           = 0.001;
  param.time           = true;
  param.test           = false;
  param.seconds        = 30;

   // for test cases

   Test *test = new Test(3,3,100000,lambda,*zX,*zY,alphX,alphY,param);
   test->comparison();
}
