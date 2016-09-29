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

 DContainer *zX = new DContainer(2,1); // alphabet
 *zX << 0 << 1;
 DContainer *zY = new DContainer(2,1); // alphabet
 *zY << 0 << 1;

 vector<double> lambda(3);
 lambda[0] = 0;
 lambda[1] = 1;
 lambda[2] = 5;

 vector<vector<int > > alphX(2,vector<int>(0));
 alphX[0].push_back(0);
 alphX[1].push_back(1);

 vector<vector<int > > alphY(2,vector<int>(0));
 alphY[0].push_back(0);
 alphY[1].push_back(1);

 vector<vector<int > > alphbX(4,vector<int>(0));
 alphbX[0].push_back(0);
 alphbX[1].push_back(0);
 alphbX[2].push_back(1);
 alphbX[3].push_back(1);

 vector<vector<int > > alphbY(4,vector<int>(0));
 alphbY[0].push_back(0);
 alphbY[1].push_back(1);
 alphbY[2].push_back(0);
 alphbY[3].push_back(1);
 cout << (*eX) << endl;

 IsParameter param;
 param.lambdavalue    = 1.0;
 param.lambdadeltaval = 1.0;
 param.sigma          = 0.01;
 param.maxit          = 10;
 param.konv           = 0.001;
 param.time           = true;
 param.test           = false;
 param.seconds        = 10;

 IT *Test = new IT(2,*eX,*zX,*zY,alphX,alphY);

 	 Test->setFeatureArraylambda(0,0,0,1);
 	 Test->setFeatureArraylambda(0,1,0,4);
 	 Test->setFeatureArraylambda(0,0,1,5);
 	 Test->setFeatureArraylambda(0,1,1,0);

 	 Test->setFeatureArraylambda(1,0,0,9);
 	 Test->setFeatureArraylambda(1,1,0,2);
 	 Test->setFeatureArraylambda(1,0,1,1);
 	 Test->setFeatureArraylambda(1,1,1,3);


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
 	 //DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, vector<double> lambda, vector<vector<int> > systAX, vector<vector<int> > systAY, vector<vector<int> > systBX, vector<vector<int> > systBY, IsParameter param, int version)
 	 TestMI *test = new TestMI(*eX,*esY,*zX,*zY,lambda,alphX, alphY,alphbX,alphbY,param,0);
 	 cout << test->getMI() << endl;
 //	 cout << (*esY) << endl;
 /*

  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;

	//V(3,vector<int>(0))
vector<vector<int > > alphX(4,vector<int>(0));
alphX[0].push_back(0);
alphX[1].push_back(1);
alphX[2].push_back(2);
alphX[3].push_back(2);

vector<vector<int > > alphY(4,vector<int>(0));
alphY[0].push_back(0);
alphY[1].push_back(1);
alphY[2].push_back(2);
alphY[3].push_back(0);
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
   test->comparison(); */
}
