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

  Test *test   = new Test(1, 1, 100, lambda, *zX, *zY, 500, 0.0001, false, true, 0, 0); // for test cases
  Test *testgp = new Test(1, 1, 100, lambda, *zX, *zY, 500, 0.0001, false, true, 0, 2);

  // p(y_0 = 0 | x_0 = 0) + p(y_0 = 1 | x_0 = 0) = 1.0
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->prop(0,   0, (*zX)(0, 0), (*zY)(0, 0))+test->prop(0,   0, (*zX)(0, 0), (*zY)(1, 0)), 0.1);
  // p(y_0 = 0 | x_0 = 1) + p(y_0 = 1 | x_0 = 1) = 1.0
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->prop(0,   0, (*zX)(1, 0), (*zY)(0, 0))+test->prop(0,   0, (*zX)(1, 0), (*zY)(1, 0)), 0.1);
  // p(y_0 = 0 | x_0 = 0) + p(y_0 = 1 | x_0 = 0) = 1.0 for test-gp
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testgp->prop(0, 0, (*zX)(0, 0), (*zY)(0, 0))+testgp->prop(0, 0, (*zX)(0, 0), (*zY)(1, 0)), 0.1);
  // p(y_0 = 0 | x_0 = 1) + p(y_0 = 1 | x_0 = 1) = 1.0 for tes-gp
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testgp->prop(0, 0, (*zX)(1, 0), (*zY)(0, 0))+testgp->prop(0, 0, (*zX)(1, 0), (*zY)(1, 0)), 0.1);

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
  CPPUNIT_ASSERT(test->KL1()   < 1.5);
  CPPUNIT_ASSERT(testgp->KL1() < 2);
}

void itTest::SCOneXOneY()
{
  cout << "SCOneXOneY" << endl;
  srand(time(NULL));
  DContainer *zX = new DContainer(2,1);
  *zX << 0 << 1;
  DContainer *zY = new DContainer(2,1);
  *zY << 0 << 1;
  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;

  Test *test   = new Test(1,1,100,lambda,*zX,*zY,500,0.0001,false,true,0,1);
  Test *testgp = new Test(1,1,100,lambda,*zX,*zY,500,0.0001,false,true,0,3);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->prop(0,   0, (*zX)(0, 0), (*zY)(0, 0))+test->prop(0,   0, (*zX)(0, 0), (*zY)(1, 0)), 0.1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->prop(0,   0, (*zX)(1, 0), (*zY)(0, 0))+test->prop(0,   0, (*zX)(1, 0), (*zY)(1, 0)), 0.1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testgp->prop(0, 0, (*zX)(0, 0), (*zY)(0, 0))+testgp->prop(0, 0, (*zX)(0, 0), (*zY)(1, 0)), 0.1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testgp->prop(0, 0, (*zX)(1, 0), (*zY)(0, 0))+testgp->prop(0, 0, (*zX)(1, 0), (*zY)(1, 0)), 0.1);

  for(int i=0; i< test->getsizeconv()-1; i++ ){
    if((test->getconv(i)>pow(10,-11))){
      if(!(test->getconv(i) >= test->getconv(i+1))){
        cout << test->getconv(i)  << " " << test->getconv(i+1) << endl;
      }
      CPPUNIT_ASSERT  (test->getconv(i) >= test->getconv(i+1) );
    }
  }
  CPPUNIT_ASSERT (test->KL1()<1.5);
  CPPUNIT_ASSERT (testgp->KL1()<2);
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

  Test *test = new Test(2,1,100,lambda,*zX,*zY,500,0.0001,false,true,0,0);
  Test *testgp = new Test(2,1,100,lambda,*zX,*zY,500,0.0001,false,true,0,2);

  for(int i=0;i<2;i++)
  {
    for(int k=0;k<2;k++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->prop(i,0,(*zX)(k,0),(*zY)(0,0))+test->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testgp->prop(i,0,(*zX)(k,0),(*zY)(0,0))+testgp->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
    }
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
  CPPUNIT_ASSERT(test->KL1()<1.5);
  CPPUNIT_ASSERT(testgp->KL1()<2);
}

void itTest::SCTwoXOneY()
{
  cout << "SCTWoXOneY" << endl;
  srand(time(NULL));
  DContainer *zX = new DContainer(2,1);
  *zX << 0 << 1;
  DContainer *zY = new DContainer(2,1);
  *zY << 0 << 1;
  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;

  Test *test   = new Test(2,1,100,lambda,*zX,*zY,500,0.0001,false,true,0,1);
  Test *testgp = new Test(2,1,100,lambda,*zX,*zY,500,0.0001,false,true,0,3);

  for(int i=0;i<2;i++)
  {
    for(int k=0;k<2;k++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->prop(i,0,(*zX)(k,0),(*zY)(0,0))+test->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testgp->prop(i,0,(*zX)(k,0),(*zY)(0,0))+testgp->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
    }
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
  CPPUNIT_ASSERT(test->KL1()<2);
  CPPUNIT_ASSERT(testgp->KL1()<2);
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

  Test *test   = new Test(2,2,100,lambda,*zX,*zY,500,0.0001,false,true,0,0);
  Test *testgp = new Test(2,2,100,lambda,*zX,*zY,500,0.0001,false,true,0,2);

  for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<2;k++){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->prop(i,0,(*zX)(k,0),(*zY)(0,0))+test->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testgp->prop(i,0,(*zX)(k,0),(*zY)(0,0))+testgp->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
      }
    }
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
  CPPUNIT_ASSERT(test->KL1()<2);
  CPPUNIT_ASSERT(testgp->KL1()<2);
}

void itTest::SCTwoXTwoY()
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

  Test *test   = new Test(2,2,100,lambda,*zX,*zY,500,0.0001,false,true,0,1);
  Test *testgp = new Test(2,2,100,lambda,*zX,*zY,500,0.0001,false,true,0,3);


  for(int i=0;i<2;i++)
  {
    for(int j=0;j<2;j++)
    {
      for(int k=0;k<2;k++)
      {
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->prop(i,0,(*zX)(k,0),(*zY)(0,0))+test->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testgp->prop(i,0,(*zX)(k,0),(*zY)(0,0))+testgp->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
      }
    }
  }

  for(int i=0; i< test->getsizeconv()-1; i++ )
  {
    if((test->getconv(i)> pow(10,-11)))
    {
      if(!(test->getconv(i) >= test->getconv(i+1)))
      {
        cout << test->getconv(i)  << " " << test->getconv(i+1) << endl;
      }
      CPPUNIT_ASSERT  (test->getconv(i) >= test->getconv(i+1) );
    }
  }
  CPPUNIT_ASSERT (test->KL1()<2);
  CPPUNIT_ASSERT (testgp->KL1()<2);
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

  Test *test = new Test(1,1,100,lambda,*zX,*zY,500,0.0001,false,true,0,0);
  Test *testgp = new Test(1,1,100,lambda,*zX,*zY,500,0.0001,false,true,0,2);

  for(int i=0;i<4;i++){
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->prop(0,0,(*zX)(i,0),(*zY)(0,0))+test->prop(0,0,(*zX)(i,0),(*zY)(1,0)),0.1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testgp->prop(0,0,(*zX)(i,0),(*zY)(0,0))+testgp->prop(0,0,(*zX)(i,0),(*zY)(1,0)),0.1);
  }


  for(int i=0; i< test->getsizeconv()-1; i++ )
  {
    CPPUNIT_ASSERT(test->getconv(i) >= test->getconv(i+1));
  }
  CPPUNIT_ASSERT (test->KL1()<1.5);
  CPPUNIT_ASSERT (testgp->KL1()<1.5);
}

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

  Test *test = new Test(4,4,100,lambda,*zX,*zY,200,0.0001,false,true,0,0);
  Test *testsc = new Test(4,4,100,lambda,*zX,*zY,200,0.0001,false,true,0,1);
  Test *testgp = new Test(4,4,100,lambda,*zX,*zY,200,0.0001,false,true,0,2);
  Test *testscgp = new Test(4,4,100,lambda,*zX,*zY,200,0.0001,false,true,0,3);

  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      for(int k=0;k<2;k++){
        CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,test->prop(i,j,k,0)+test->prop(i,j,k,1)+test->prop(i,j,k,2),0.1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,testsc->prop(i,j,k,0)+testsc->prop(i,j,k,1)+testsc->prop(i,j,k,2),0.1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,testgp->prop(i,j,k,0)+testgp->prop(i,j,k,1)+testgp->prop(i,j,k,2),0.1);
        CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,testscgp->prop(i,j,k,0)+testscgp->prop(i,j,k,1)+testscgp->prop(i,j,k,2),0.1);
      }
    }
  }
  CPPUNIT_ASSERT (test->KL1()<2);
  CPPUNIT_ASSERT (testsc->KL1()<2);
  CPPUNIT_ASSERT (testgp->KL1()<2);
  CPPUNIT_ASSERT (testscgp->KL1()<2);
}
