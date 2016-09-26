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
  vector<int> x(3);
  x[0]=0;
  x[1]=1;
  x[2]=2;
  vector<int> y(2);
  y[0]=0;
  y[1]=1;
  //DContainer &aX, DContainer &aY,vector<double> systX,vector<double> systY, double valuelambda
 // Feature *M = new Feature(*zX,*zY,3,2,x,y,2);
 // cout << (*M)<< endl;

  //int indexX, int indexY, vector<double> x, vector<double> y
  vector<double> xi(3);
  xi[0]=0;
  xi[1]=0;
  xi[2]=0;
  cout << " x3 " << xi[2] << endl;
  vector<double> yj(2);
  yj[0]=0;
  yj[1]=0;

  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.001;
  param.time           = true;
  param.test           = false;
  param.seconds        = 20;

//  Test *test = new Test(2, 2, 100000, *zX, *zY, param); // for test cases
  // Test *test = new Test(5, 5, 100000, lambda, *zX, *zY, 500, 0.0001, true, 2);
//  test->comparison();
  	 srand(time(NULL));
  int n=100;
  DContainer *eX = new DContainer(n,2);
  for(int i=0;i< n;i++ ){
  	for(int j=0;j<2;j++){
  		*eX << rand() % 2;
 	}
  }
 // DContainer *eX = new DContainer(10,2);
 // (*eX) << 0 << 0 << 0 << 0 << 0 <<1 << 1 << 1 << 1 << 1  << 1 << 1 << 1 << 1 << 1 ;
  cout << (*eX) << endl;

  	//V(3,vector<int>(0))
 vector<vector<int > > alphX(4,vector<int>(0));
 alphX[0].push_back(0);
 alphX[1].push_back(0);
 alphX[2].push_back(1);
 alphX[3].push_back(1);

 vector<vector<int > > alphY(4,vector<int>(0));
 alphY[0].push_back(0);
 alphY[1].push_back(1);
 alphY[2].push_back(0);
 alphY[3].push_back(1);

  IT *Test = new IT(2,*eX,*zX,*zY,alphX,alphY);

  	 Test->setFeatureArraylambda(0,0,0,1);
  	 Test->setFeatureArraylambda(0,1,0,4);
  	 Test->setFeatureArraylambda(0,0,1,5);
  	 Test->setFeatureArraylambda(0,1,1,0);

  	 Test->setFeatureArraylambda(1,0,0,0);
  	 Test->setFeatureArraylambda(1,1,0,2);
  	 Test->setFeatureArraylambda(1,0,1,1);
  	 Test->setFeatureArraylambda(1,1,1,3);

  	 Test->setFeatureArraylambda(2,0,0,3);
  	 Test->setFeatureArraylambda(2,1,0,2);
  	 Test->setFeatureArraylambda(2,0,1,0);
  	 Test->setFeatureArraylambda(2,1,1,3);

  	 Test->setFeatureArraylambda(3,0,0,2);
  	 Test->setFeatureArraylambda(3,1,0,1);
  	 Test->setFeatureArraylambda(3,0,1,0);
  	 Test->setFeatureArraylambda(3,1,1,3);

  	 double** prop;
  	 prop=new double*[n];
  	 for(int i=0;i<n;i++){
  		 prop[i]=new double[4];
  		 for(int j=0;j<4;j++){
  			 prop[i][j]=0;
  		 }
  	 }

  	for(int i=0;i<n;i++){
  		for(int propi=0;propi<4;propi++){
  				prop[i][propi]=Test->prop(i,propi);
  			}
  	}
  	 DContainer *esY=new DContainer(n,2);

  	 for(int i=0;i<n;i++){
  		 double z=(double)rand()/RAND_MAX;
  		 double s=0;
  		 int ind=0;
  		 for(int j=0;j<4 && s<z;j++){
  				s+=prop[i][j];
  				ind=j;
  		 }
  		 if(ind==0) (*esY) << 0 << 0;
  		 if(ind==1) (*esY) << 1 << 0;
  		 if(ind==2) (*esY) << 0 << 1;
  		 if(ind==3) (*esY) << 1 << 1;

  	 }
  	 cout << (*esY) << endl;
  	// cout << (*esY) << endl;
  	 //DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param
 	 SCGIS *_scgisTest = new SCGIS(*eX,*esY,*zX,*zY,alphX,alphY,param);
  	 GIS *_gisTest = new GIS(*eX,*esY,*zX,*zY,alphX,alphY,param);
     cout << "lambda: " << endl;
     cout << "vergleichswerte" << endl;
     cout << Test->getFeatureArraylambda(0,0,0) <<endl;
     cout << Test->getFeatureArraylambda(0,1,0) <<endl;
     cout << Test->getFeatureArraylambda(0,0,1) <<endl;
     cout << Test->getFeatureArraylambda(0,1,1) <<endl;
     cout << endl;

     cout << Test->getFeatureArraylambda(1,0,0) <<endl;
     cout << Test->getFeatureArraylambda(1,1,0) <<endl;
     cout << Test->getFeatureArraylambda(1,0,1) <<endl;
     cout << Test->getFeatureArraylambda(1,1,1) <<endl;
     cout << endl;
     cout << Test->getFeatureArraylambda(2,0,0) <<endl;
     cout << Test->getFeatureArraylambda(2,1,0) <<endl;
     cout << Test->getFeatureArraylambda(2,0,1) <<endl;
     cout << Test->getFeatureArraylambda(2,1,1) <<endl;
     cout << endl;
     cout << Test->getFeatureArraylambda(3,0,0) <<endl;
     cout << Test->getFeatureArraylambda(3,1,0) <<endl;
     cout << Test->getFeatureArraylambda(3,0,1) <<endl;
     cout << Test->getFeatureArraylambda(3,1,1) <<endl;

     cout << "GIS" << endl;
     cout << _gisTest->getFeatureArraylambda(0,0,0) <<endl;
     cout << _gisTest->getFeatureArraylambda(0,1,0) <<endl;
     cout << _gisTest->getFeatureArraylambda(0,0,1) <<endl;
     cout << _gisTest->getFeatureArraylambda(0,1,1) <<endl;
     cout << endl;

     cout << _gisTest->getFeatureArraylambda(1,0,0) <<endl;
     cout << _gisTest->getFeatureArraylambda(1,1,0) <<endl;
     cout << _gisTest->getFeatureArraylambda(1,0,1) <<endl;
     cout << _gisTest->getFeatureArraylambda(1,1,1) <<endl;
     cout << endl;
     cout << _gisTest->getFeatureArraylambda(2,0,0) <<endl;
     cout << _gisTest->getFeatureArraylambda(2,1,0) <<endl;
     cout << _gisTest->getFeatureArraylambda(2,0,1) <<endl;
     cout << _gisTest->getFeatureArraylambda(2,1,1) <<endl;
     cout << endl;
     cout << _gisTest->getFeatureArraylambda(3,0,0) <<endl;
     cout << _gisTest->getFeatureArraylambda(3,1,0) <<endl;
     cout << _gisTest->getFeatureArraylambda(3,0,1) <<endl;
     cout << _gisTest->getFeatureArraylambda(3,1,1) <<endl;
     cout << "SCGIS" << endl;
    cout << _scgisTest->getFeatureArraylambda(0,0,0) <<endl;
     cout << _scgisTest->getFeatureArraylambda(0,1,0) <<endl;
     cout << _scgisTest->getFeatureArraylambda(0,0,1) <<endl;
     cout << _scgisTest->getFeatureArraylambda(0,1,1) <<endl;
     cout << endl;

     cout << _scgisTest->getFeatureArraylambda(1,0,0) <<endl;
     cout << _scgisTest->getFeatureArraylambda(1,1,0) <<endl;
     cout << _scgisTest->getFeatureArraylambda(1,0,1) <<endl;
     cout << _scgisTest->getFeatureArraylambda(1,1,1) <<endl;
     cout << endl;
     cout << _scgisTest->getFeatureArraylambda(2,0,0) <<endl;
     cout << _scgisTest->getFeatureArraylambda(2,1,0) <<endl;
     cout << _scgisTest->getFeatureArraylambda(2,0,1) <<endl;
     cout << _scgisTest->getFeatureArraylambda(2,1,1) <<endl;
     cout << endl;
     cout << _scgisTest->getFeatureArraylambda(3,0,0) <<endl;
     cout << _scgisTest->getFeatureArraylambda(3,1,0) <<endl;
     cout << _scgisTest->getFeatureArraylambda(3,0,1) <<endl;
     cout << _scgisTest->getFeatureArraylambda(3,1,1) <<endl;
}
