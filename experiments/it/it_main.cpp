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
  DContainer *zX = new DContainer(3,1); // alphabet
  *zX << 0 << 1 << 2;
  DContainer *zY = new DContainer(2,1); // alphabet
  *zY << 0 << 1;
  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;
  vector<int> x(3);
  x[0]=0;
  x[1]=2;
  x[2]=3;
  vector<int> y(2);
  y[0]=0;
  y[1]=1;
  //DContainer &aX, DContainer &aY,vector<double> systX,vector<double> systY, double valuelambda
  Feature *M = new Feature(*zX,*zY,x,y,2);
  cout << (*M)<< endl;
  for(int i=0;i<27;i++){
	  vector<double> ind = M->index(i,true);
  }
  //int indexX, int indexY, vector<double> x, vector<double> y
  vector<double> xi(3);
  xi[0]=0;
  xi[1]=0;
  xi[2]=1;
 // cout << " x3 " << x[2] << endl;
  vector<double> yj(2);
  yj[0]=0;
  yj[1]=0;
  cout << " delta " << M->delta(1,0,xi,yj) << endl;
  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.001;
  param.time           = true;
  param.test           = false;
  param.seconds        = 30;

//  Test *test = new Test(2, 2, 100000, *zX, *zY, param); // for test cases
  // Test *test = new Test(5, 5, 100000, lambda, *zX, *zY, 500, 0.0001, true, 2);
//  test->comparison();
}
