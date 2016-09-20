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



int main(int argc, char **argv){

DContainer *zX = new DContainer(2,1);
  *zX << 0 << 1;
DContainer *zY = new DContainer(2,1);
  *zY << 0 << 1;

  vector<double> lambda(3);
     lambda[0]=1;
     lambda[1]=5;
     lambda[2]=3;
//int ColX,int RowX,int ColValY,  vector<double> lambda,DContainer &aX, DContainer &aY,int maxit, double konv,bool time,int seconds)
     Test *test = new Test(2,2,100,lambda,*zX,*zY,1,0.001,true,10);
test->comparison();
}
