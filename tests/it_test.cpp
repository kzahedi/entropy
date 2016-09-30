#include "it_test.h"

#include <entropy++/Container.h>
#include <entropy++/MI.h>
#include <entropy++/sparse/MI.h>
#include "../experiments/it/GIS.h"
#include "../experiments/it/SCGIS.h"
#include "../experiments/it/Test.h"

#include <iostream>
#include <string>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( itTest );

void itTest::OneXOneY()
{
  cout << "OneXOneY" << endl;
  srand(time(NULL));
  DContainer *zX = new DContainer(2,1); // alphabet
  *zX << 0 << 1;
  DContainer *zY = new DContainer(2,1); // alphabet
  *zY << 0 << 1;
  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;

  vector<vector<int > > alphX(1,vector<int>(0));
  alphX[0].push_back(0);

  vector<vector<int > > alphY(1,vector<int>(0));
  alphY[0].push_back(0);

  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.0001;
  param.time           = false;
  param.test           = true;
  param.seconds        = 0;

  IsParameter paramgp;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.0001;
  param.time           = true;
  param.test           = true;
  param.seconds        = 5;

  Test *testval = new Test(1,1,100,lambda, *zX,*zY,alphX, alphY); //get data

  GIS *test         = new GIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  GISgp *testgp     = new GISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  SCGISgp *testscgp = new SCGISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);

  // p(y_0 = 0 | x_0 = 0) + p(y_0 = 1 | x_0 = 0) = 1.0
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->prop( 0, 0)+test->prop( 0, 1), 0.1);
  // p(y_0 = 0 | x_0 = 1) + p(y_0 = 1 | x_0 = 1) = 1.0
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->prop( 1, 0)+test->prop( 1, 1), 0.1);
  // p(y_0 = 0 | x_0 = 0) + p(y_0 = 1 | x_0 = 0) = 1.0 for test-gp
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testgp->prop(0, 0)+testgp->prop(0,1), 0.1);
  // p(y_0 = 0 | x_0 = 1) + p(y_0 = 1 | x_0 = 1) = 1.0 for tes-gp
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testgp->prop(0, 0)+testgp->prop(0,1), 0.1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testsc->prop(0, 0)+testsc->prop(0, 1), 0.1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testsc->prop(1, 0)+testsc->prop(1, 1), 0.1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testscgp->prop(0, 0)+testscgp->prop(0, 1), 0.1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testscgp->prop(1, 0)+testscgp->prop(1, 1), 0.1);

  // getsizeconv -> size of stored l-values (constraints)
  // constraints: difference betweens observed feature frequency vs. expected feature frequency, given the lambdas
  for(int i = 0; i < test->getsizeconv()-1; i++)
  {
    if((test->getconv(i)>pow(10,-11)))
    {
      if(!(test->getconv(i) >= test->getconv(i+1))) // must be decreasing, for output
      {
        cout << test->getconv(i)  << " " << test->getconv(i+1) << endl;
      }
      CPPUNIT_ASSERT  (test->getconv(i) >= test->getconv(i+1) );
    }
  }
  for(int i=0; i< testsc->getsizeconv()-1; i++ ){
      if((testsc->getconv(i)>pow(10,-11))){
        if(!(testsc->getconv(i) >= testsc->getconv(i+1))){
          cout << testsc->getconv(i)  << " " << testsc->getconv(i+1) << endl;
        }
        CPPUNIT_ASSERT  (testsc->getconv(i) >= testsc->getconv(i+1) );
      }
    }
 // cout <<"KL1 " <<  test->KL1() << endl;
 // cout <<"KL1gp " <<  testgp->KL1() << endl;
//  CPPUNIT_ASSERT(test->KL1()   < 1);
//  CPPUNIT_ASSERT(testgp->KL1() < 1);
}
void itTest::TwoXOneY()
{
  cout << "TWoXOneY" << endl;
  srand(time(NULL));
  DContainer *zX = new DContainer(2,1);
  *zX << 0 << 1;
  DContainer *zY = new DContainer(2,1);
  *zY << 0 << 1;
  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;

  vector<vector<int > > alphX(3,vector<int>(0));
  alphX[0].push_back(0);
  alphX[0].push_back(1);
  alphX[1].push_back(0);
  alphX[2].push_back(1);

  vector<vector<int > > alphY(3,vector<int>(0));
  alphY[0].push_back(0);
  alphY[1].push_back(0);
  alphY[2].push_back(0);

  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.0001;
  param.time           = false;
  param.test           = true;
  param.seconds        = 0;

  IsParameter paramgp;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.0001;
  param.time           = true;
  param.test           = true;
  param.seconds        = 5;

  Test *testval   = new Test(2,1,100,lambda, *zX,*zY,alphX, alphY);
  GIS *test     = new GIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  GISgp *testgp = new GISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  SCGISgp *testscgp = new SCGISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);

    for(int k=0;k<2;k++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->prop(k,0)+test->prop(k,1),0.1);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testgp->prop(k,0)+testgp->prop(k,1),0.1);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testsc->prop(k,0)+testsc->prop(k,1),0.1);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testscgp->prop(k,0)+testscgp->prop(k,1),0.1);
    }

  for(int i=0; i< test->getsizeconv()-1; i++ )
  {
    if((test->getconv(i)>pow(10,-11)))
    {
      if(!(test->getconv(i) >= test->getconv(i+1)))
      {
        cout << test->getconv(i)  << " " << test->getconv(i+1) << endl;
      }
      CPPUNIT_ASSERT(test->getconv(i) >= test->getconv(i+1));
    }
  }
  for(int i=0; i< testsc->getsizeconv()-1; i++ )
  {
    if((testsc->getconv(i)>pow(10,-11)))
    {
      if(!(testsc->getconv(i) >= testsc->getconv(i+1)))
      {
        cout << testsc->getconv(i)  << " " << testsc->getconv(i+1) << endl;
      }
      CPPUNIT_ASSERT  (testsc->getconv(i) >= testsc->getconv(i+1) );
    }
  }
 // cout <<"KL1 " <<  test->KL1() << endl;
 // cout << "KL1gp" << testgp->KL1() << endl;
  //CPPUNIT_ASSERT(test->KL1()<1);
//  CPPUNIT_ASSERT(testgp->KL1()<2);
}
void itTest::TwoXTwoY()
{
  cout << "TwoXTwoY" << endl;
  srand(time(NULL));
  DContainer *zX = new DContainer(2,1);
  *zX << 0 << 1;
  DContainer *zY = new DContainer(2,1);
  *zY << 0 << 1;
  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;

  vector<vector<int > > alphX(3,vector<int>(0));
  alphX[0].push_back(0);
  alphX[0].push_back(1);
  alphX[1].push_back(0);
  alphX[2].push_back(1);

  vector<vector<int > > alphY(3,vector<int>(0));
  alphY[0].push_back(0);
  alphY[1].push_back(1);
  alphY[2].push_back(0);

  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.0001;
  param.time           = false;
  param.test           = true;
  param.seconds        = 0;

  IsParameter paramgp;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.0001;
  param.time           = true;
  param.test           = true;
  param.seconds        = 5;

  Test *testval   = new Test(2,2,100,lambda, *zX,*zY,alphX, alphY);
  GIS *test     = new GIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  GISgp *testgp = new GISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  SCGISgp *testscgp = new SCGISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);

  for(int i=0;i<4;i++){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->prop(i,0)+test->prop(i,1)+ test->prop(i,2) +test->prop(i,3),0.1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testgp->prop(i,0)+testgp->prop(i,1)+ testgp->prop(i,2) + testgp->prop(i,3),0.1);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testsc->prop(i,0)+test->prop(i,1)+ testsc->prop(i,2) +test->prop(i,3),0.1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testscgp->prop(i,0)+testgp->prop(i,1)+ testscgp->prop(i,2) +testgp->prop(i,3),0.1);
  }

  for(int i=0; i< test->getsizeconv()-1; i++ )
  {
    if((test->getconv(i)>pow(10,-11)))
    {
      if(!(test->getconv(i) >= test->getconv(i+1)))
      {
        cout << test->getconv(i)  << " " << test->getconv(i+1) << endl;
      }
      CPPUNIT_ASSERT  (test->getconv(i) >= test->getconv(i+1) );
    }
  }

  for(int i=0; i< testsc->getsizeconv()-1; i++ )
  {
    if((testsc->getconv(i)> pow(10,-11)))
    {
      if(!(testsc->getconv(i) >= testsc->getconv(i+1)))
      {
        cout << testsc->getconv(i)  << " " << testsc->getconv(i+1) << endl;
      }
      CPPUNIT_ASSERT  (testsc->getconv(i) >= testsc->getconv(i+1) );
    }
  }
 // CPPUNIT_ASSERT(test->KL1()<1);
 // CPPUNIT_ASSERT(testgp->KL1()<2);
//  cout << "KL1 "<< test->KL1() << endl;
 // cout << "KL1gp " <<testgp->KL1() << endl;
//  cout << testgp->KL1() << endl;
}

void itTest::NotBinary()
{
  cout << "NotBinary" << endl;
  srand(time(NULL));
  DContainer *zX = new DContainer(4,1);
  *zX << 0 << 1 << 2 << 3;
  DContainer *zY = new DContainer(2,1);
  *zY << 0 << 1;
  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;

  vector<vector<int > > alphX(1,vector<int>(0));
  alphX[0].push_back(0);


  vector<vector<int > > alphY(1,vector<int>(0));
  alphY[0].push_back(0);

  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.0001;
  param.time           = false;
  param.test           = true;
  param.seconds        = 0;

  IsParameter paramgp;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.0001;
  param.time           = true;
  param.test           = true;
  param.seconds        = 5;

  Test *testval   = new Test(1,1,100,lambda, *zX,*zY,alphX, alphY);
  GIS *test     = new GIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  GISgp *testgp = new GISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);

  for(int i=0;i<4;i++){
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->prop(i,0)+test->prop(i,1),0.1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testgp->prop(i,0)+testgp->prop(i,1),0.1);
  }


  for(int i=0; i< test->getsizeconv()-1; i++ )
  {
    CPPUNIT_ASSERT(test->getconv(i) >= test->getconv(i+1));
  }
//  cout <<"KL1 " <<  test->KL1() << endl;
//  cout <<"KL1gp " << testgp->KL1() << endl;
 // CPPUNIT_ASSERT (test->KL1()<1);
 // CPPUNIT_ASSERT (testgp->KL1()<1.5);
}
/*
void itTest::FourXFourY()
{
  cout << "FourXFourY" << endl;
  srand(time(NULL));
  DContainer *zX = new DContainer(2,1);
  *zX << 0 << 1;
  DContainer *zY = new DContainer(3,1);
  *zY << 0 << 1 << 2;
  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;

  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 1;
  param.konv           = 0.0001;
  param.time           = true;
  param.test           = true;
  param.seconds        = 5;


  vector<vector<int > > alphX(4,vector<int>(0));
  alphX[0].push_back(0);
  alphX[0].push_back(1);
  alphX[1].push_back(0);
  alphX[2].push_back(1);
  alphX[2].push_back(3);
  alphX[3].push_back(2);

  vector<vector<int > > alphY(4,vector<int>(0));
  alphY[0].push_back(0);
  alphY[1].push_back(1);
  alphY[2].push_back(0);
  alphY[3].push_back(2);
//ohne cases...
  Test *test     = new Test(4,4,100,lambda,*zX,*zY, alphX, alphY, param);

  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  double sum4 = 0.0;

      for(int k=0;k<16;k++){
    	  for(int i=0;i<81;i++){
   		  sum1+=test->propAll(k,i)[0];
   		  sum2+=test->propAll(k,i)[1];
    		  sum3+=test->propAll(k,i)[2];
   		  sum4+=test->propAll(k,i)[3];
    	  }
          CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum1,0.1);
          CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum2,0.1);
          CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum3,0.1);
          CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum4,0.1);
          sum1 = 0;
          sum2 = 0;
          sum3 = 0;
          sum4 = 0;
  }
      cout <<"KL1 " <<  test->KL()[0] << endl;
      cout <<"KL1sc " <<  test->KL()[1]<< endl;
      cout <<"KL1gp"  << test->KL()[2] << endl;
      cout <<"KL1scgp " << test->KL()[3]<< endl;
 // CPPUNIT_ASSERT (test->KL1()<2);
 // CPPUNIT_ASSERT (testsc->KL1()<2);
 // CPPUNIT_ASSERT (testgp->KL1()<2);
 // CPPUNIT_ASSERT (testscgp->KL1()<2);
}

*/
