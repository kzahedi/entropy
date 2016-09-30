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
#include "TestMI.h"



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

 vector<vector<int > > alphbX(4,vector<int>(0));
 alphbX[0].push_back(0);
 alphbX[0].push_back(1);

 vector<vector<int > > alphbY(4,vector<int>(0));
 alphbY[0].push_back(0);


 IsParameter param;
 param.lambdavalue    = 1.0;
 param.lambdadeltaval = 1.0;
 param.sigma          = 0.01;   //TODO find best sigma value
 param.maxit          = 10;
 param.konv           = 0.000001;
 param.time           = true;
 param.test           = true;
 param.seconds        = 10;

 vector<int> cases(4);
 cases[0]=0;
 cases[1]=1;
 cases[2]=2;
 cases[3]=3;
  Test *test = new Test(2,1,100,lambda, *zX,*zY,alphbX, alphbY);
  test->compareCases(param,cases);
  //(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, IsParameter param, int version)
  TestMI *testmi = new TestMI(test->getvalX(),test->getvalY(),*zX,*zY,param,0);
  cout << testmi->getMI() << endl;
  TestMI *testmisc = new TestMI(test->getvalX(),test->getvalY(),*zX,*zY,param,1);
  cout << testmisc->getMI() << endl;
  TestMI *testmigp = new TestMI(test->getvalX(),test->getvalY(),*zX,*zY,param,2);
  cout << testmigp->getMI() << endl;
  TestMI *testmiscgp = new TestMI(test->getvalX(),test->getvalY(),*zX,*zY,param,3);
  cout << testmiscgp->getMI() << endl;
}
