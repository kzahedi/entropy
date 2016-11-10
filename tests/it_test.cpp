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
void itTest::ValueTest()
{  DContainer *zX = new DContainer(2,1); // alphabet
      *zX << 0 << 1;
   DContainer *zY = new DContainer(2,1); // alphabet
      *zY << 0 << 1;
   vector<vector<int > > alphX(1,vector<int>(0));
      alphX[0].push_back(0);
   vector<vector<int > > alphY(1,vector<int>(0));
      alphY[0].push_back(0);
   vector<double> lambda(1);
       lambda[0] = 1;
   Test *test_1 = new Test(1,1,100,lambda, *zX,*zY,alphX, alphY);
   for(int i=0; i<2;i++){
	   for(int j=0; j<2;j++){
          CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, test_1->getProp(i,j), 0.1);
	   }
   }
   IContainer *indizes = new IContainer(1,3);
   (*indizes) << 0 << 0 << 0;
   DContainer *values  = new DContainer(1,1);
   (*values) << 5;
   Test *test_2 = new Test(1,1,100,*indizes,*values, *zX,*zY,alphX, alphY);
   CPPUNIT_ASSERT( test_2->getProp(0,0)> test_2->getProp(0,1));
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, test_2->getProp(1,1),0.01);
   CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, test_2->getProp(1,0),0.01);

}
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
  param.konv           = 0.00001;
  param.konvtime       = true;
  param.time           = false;
  param.test           = true;
  param.seconds        = 72000;
/*
  IsParameter paramgp;
  paramgp.lambdavalue    = 1.0;
  paramgp.lambdadeltaval = 1.0;
  paramgp.sigma          = 0.01;
  paramgp.maxit          = 500;
  paramgp.konv           = 0.0001;
  paramgp.time           = true;
  paramgp.test           = true;
  paramgp.seconds        = 60; */

  Test *testval = new Test(1,1,10000,lambda, *zX,*zY,alphX, alphY); //get data

  GIS *test         = new GIS(testval->getvalX(),testval->getvalY(),*zX,*zY,alphX,alphY,param);
//  GISgp *testgp     = new GISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,paramgp);
  SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zX,*zY,alphX,alphY,param);
//  SCGISgp *testscgp = new SCGISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,paramgp);

  // p(y_0 = 0 | x_0 = 0) + p(y_0 = 1 | x_0 = 0) = 1.0
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->propAlphX( 0, 0)+test->propAlphX( 0, 1), 0.1);
  // p(y_0 = 0 | x_0 = 1) + p(y_0 = 1 | x_0 = 1) = 1.0
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->propAlphX( 1, 0)+test->propAlphX( 1, 1), 0.1);
  // p(y_0 = 0 | x_0 = 0) + p(y_0 = 1 | x_0 = 0) = 1.0 for test-gp
//  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testgp->propAlphX(0, 0)+testgp->propAlphX(0,1), 0.1);
  // p(y_0 = 0 | x_0 = 1) + p(y_0 = 1 | x_0 = 1) = 1.0 for tes-gp
//  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testgp->propAlphX(0, 0)+testgp->propAlphX(0,1), 0.1);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testsc->propAlphX(0, 0)+testsc->propAlphX(0, 1), 0.1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testsc->propAlphX(1, 0)+testsc->propAlphX(1, 1), 0.1);
//  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testscgp->propAlphX(0, 0)+testscgp->propAlphX(0, 1), 0.1);
//  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testscgp->propAlphX(1, 0)+testscgp->propAlphX(1, 1), 0.1);

  double exact;
  for(int i=0;i< 2; i++){
	  for(int j=0; j< 2; j++){
		  exact= testval->getProp(i,j);
//		  CPPUNIT_ASSERT_DOUBLES_EQUAL(exact,test->propAlphX(i,j),0.2);
//		  CPPUNIT_ASSERT_DOUBLES_EQUAL(exact,testsc->propAlphX(i,j),0.2);
	//	  CPPUNIT_ASSERT_DOUBLES_EQUAL(exact,testgp->propAlphX(i,j),0.2);
	//	  CPPUNIT_ASSERT_DOUBLES_EQUAL(exact,testscgp->propAlphX(i,j),0.2);
		  cout << i << " " << j << endl;
		  cout << "exakt: " <<  testval->getProp(i,j) << endl;
		  cout << "gis:  " << test->propAlphX(i,j) << endl;
		  cout << "scgis: " << testsc->propAlphX(i,j) << endl;
	//	  cout << "gisgp: " << testgp->propAlphX(i,j) << endl;
	//	  cout << "scgisgp: " << testscgp->propAlphX(i,j) << endl;
	  }
  }



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
  cout << " GIS: " << testval->KL(test);
  cout << " SCGIS: " << testval->KL(testsc);
//  cout << " GISgp: " << testval->KL(testgp);
//  cout << " SCGISgp: " << testval->KL(testscgp);
}
void itTest:: OneXOneYFour(){
	  cout << "OneXOneYFour" << endl;
	  srand(time(NULL));
	  DContainer *zX = new DContainer(4,1); // alphabet
	  *zX << 0 << 1 << 2 << 3;
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
	  param.konvtime       = true;
	  param.time           = false;
	  param.test           = true;
	  param.seconds        = 72000;


	  Test *testval = new Test(1,1,10000,lambda, *zX,*zX,alphX, alphY); //get data

	  GIS *test         = new GIS(testval->getvalX(),testval->getvalY(),*zX,*zX,alphX,alphY,param);
	  SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zX,*zX,alphX,alphY,param);

	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->propAlphX( 0, 0)+test->propAlphX( 0, 1)+test->propAlphX( 0, 2)+test->propAlphX( 0, 3), 0.1);
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->propAlphX( 1, 0)+test->propAlphX( 1, 1)+test->propAlphX( 1, 2)+test->propAlphX( 1, 3), 0.1);
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->propAlphX( 2, 0)+test->propAlphX( 2, 1)+test->propAlphX( 2, 2)+test->propAlphX( 2, 3), 0.1);
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, test->propAlphX( 3, 0)+test->propAlphX( 3, 1)+test->propAlphX( 3, 2)+test->propAlphX( 3, 3), 0.1);

	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testsc->propAlphX( 0, 0)+testsc->propAlphX( 0, 1)+testsc->propAlphX( 0, 2)+testsc->propAlphX( 0, 3), 0.1);
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testsc->propAlphX( 1, 0)+testsc->propAlphX( 1, 1)+testsc->propAlphX( 1, 2)+testsc->propAlphX( 1, 3), 0.1);
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testsc->propAlphX( 2, 0)+testsc->propAlphX( 2, 1)+testsc->propAlphX( 2, 2)+testsc->propAlphX( 2, 3), 0.1);
	  CPPUNIT_ASSERT_DOUBLES_EQUAL(1, testsc->propAlphX( 3, 0)+testsc->propAlphX( 3, 1)+testsc->propAlphX( 3, 2)+testsc->propAlphX( 3, 3), 0.1);

	  double exact;
	  for(int i=0;i< 4; i++){
		  for(int j=0; j< 4; j++){
			  cout << i << " " << j << endl;
			  cout << "exakt: " <<  testval->getProp(i,j) << endl;
			  cout << "gis:  " << test->propAlphX(i,j) << endl;
			  cout << "scgis: " << testsc->propAlphX(i,j) << endl;
		  }
	  }
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
	  cout << " GIS: " << testval->KL(test);
	  cout << " SCGIS: " << testval->KL(testsc);
}
void itTest::OneXOneYEight(){
	  cout << "OneXOneYEight" << endl;
	  srand(time(NULL));
	  DContainer *zX = new DContainer(8,1); // alphabet
	  *zX << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7;
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
	  param.konvtime       = true;
	  param.time           = false;
	  param.test           = true;
	  param.seconds        = 72000;


	  Test *testval = new Test(1,1,10000,lambda, *zX,*zX,alphX, alphY); //get data

	  GIS *test         = new GIS(testval->getvalX(),testval->getvalY(),*zX,*zX,alphX,alphY,param);
	  SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zX,*zX,alphX,alphY,param);

	  double sum1 = 0.0;
	  double sum2 = 0.0;

	      for(int k=0;k<8;k++){
	    	  for(int i=0;i<8;i++){
	   		    sum1+=test->propAlphX(k,i);
	            sum2+=testsc->propAlphX(k,i);
	    	  }
	          CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum1,0.1);
	          CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum2,0.1);
	          sum1 = 0;
	          sum2 = 0;
	  }

	  double exact;
	  for(int i=0;i< 8; i++){
		  for(int j=0; j< 8; j++){
			  cout << i << " " << j << endl;
			  cout << "exakt: " <<  testval->getProp(i,j) << endl;
			  cout << "gis:  " << test->propAlphX(i,j) << endl;
			  cout << "scgis: " << testsc->propAlphX(i,j) << endl;
		  }
	  }
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
	  cout << " GIS: " << testval->KL(test);
	  cout << " SCGIS: " << testval->KL(testsc);
}
void itTest::OneXOneYTwenty(){
	  cout << "OneXOneYTwenty" << endl;
	  srand(time(NULL));
	  DContainer *zX = new DContainer(20,1); // alphabet
	  for(int i=0;i<20;i++){
		  (*zX) <<i;
	  }
	  vector<double> lambda(3);
	  lambda[0] = 0;
	  lambda[1] = 1;
	  lambda[2] = 2;

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
	  param.konvtime       = true;
	  param.time           = false;
	  param.test           = true;
	  param.seconds        = 72000;


	  Test *testval = new Test(1,1,10000,lambda, *zX, *zX, alphX, alphY); //get data

	  GIS *test         = new GIS(testval->getvalX(),testval->getvalY(),*zX,*zX,alphX,alphY,param);
	  SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zX,*zX,alphX,alphY,param);

	  double sum1 = 0.0;
	  double sum2 = 0.0;

	      for(int k=0;k<20;k++){
	    	  for(int i=0;i<20;i++){
	   		    sum1+=test->propAlphX(k,i);
	            sum2+=testsc->propAlphX(k,i);
	    	  }
	          CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum1,0.1);
	          CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum2,0.1);
	          sum1 = 0;
	          sum2 = 0;
	  }

	  for(int i=0;i<20; i++){
		  for(int j=0; j<20; j++){
			  cout << i << " " << j << endl;
			  cout << "exakt: " <<  testval->getProp(i,j) << endl;
			  cout << "gis:  " << test->propAlphX(i,j) << endl;
			  cout << "scgis: " << testsc->propAlphX(i,j) << endl;
		  }
	  }
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
	  cout << " GIS: " << testval->KL(test);
	  cout << " SCGIS: " << testval->KL(testsc);
}
void itTest::TwoXOneY()
{
  cout << endl;
  cout << "TWoXOneY" << endl;
  srand(time(NULL));
  DContainer *zY = new DContainer(2,1);
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
  alphY[1].push_back(0);

  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.0001;
  param.konvtime       = true;
  param.time           = false;
  param.test           = true;
  param.seconds        = 72000;

  /*
  IsParameter paramgp;
  paramgp.lambdavalue    = 1.0;
  paramgp.lambdadeltaval = 1.0;
  paramgp.sigma          = 0.01;
  paramgp.maxit          = 500;
  paramgp.konv           = 0.0001;
  paramgp.time           = true;
  paramgp.test           = true;
  paramgp.seconds        = 60; */

  Test *testval     = new Test(2,1,10000,lambda, *zY,*zY,alphX, alphY);
  GIS *test         = new GIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
//  GISgp *testgp     = new GISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,paramgp);
  SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
//  SCGISgp *testscgp = new SCGISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,paramgp);

    for(int k=0;k<2;k++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->propAlphX(k,0)+test->propAlphX(k,1),0.1);
//      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testgp->propAlphX(k,0)+testgp->propAlphX(k,1),0.1);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testsc->propAlphX(k,0)+testsc->propAlphX(k,1),0.1);
//      CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testscgp->propAlphX(k,0)+testscgp->propAlphX(k,1),0.1);
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
  double exact;
  for(int i=0;i< 4; i++){
	  for(int j=0; j< 2; j++){
		  cout << i << " " << j << endl;
		  cout << "exakt: " <<  testval->getProp(i,j) << endl;
		  cout << "gis:  " << test->propAlphX(i,j) << endl;
		  cout << "scgis: " << testsc->propAlphX(i,j) << endl;
	//	  cout << "gisgp: " << testgp->propAlphX(i,j) << endl;
	//	  cout << "scgisgp: " << testscgp->propAlphX(i,j) << endl;
		  exact= testval->getProp(i,j);
    //	  CPPUNIT_ASSERT_DOUBLES_EQUAL(exact,test->propAlphX(i,j),0.3);
    //    CPPUNIT_ASSERT_DOUBLES_EQUAL(exact,testsc->propAlphX(i,j),0.3);
	//	  CPPUNIT_ASSERT_DOUBLES_EQUAL(exact,testgp->propAlphX(i,j),0.3);
	//	  CPPUNIT_ASSERT_DOUBLES_EQUAL(exact,testscgp->propAlphX(i,j),0.3);
	  }
  }

  cout << " GIS: " << testval->KL(test);
  cout << " SCGIS: " << testval->KL(testsc);
//  cout << " GISgp: " << testval->KL(testgp);
//  cout << " SCGISgp: " << testval->KL(testscgp);
}
void itTest::TwoXTwoY()
{
  cout << endl;
  cout << "TwoXTwoY" << endl;
  srand(time(NULL));
  DContainer *zY = new DContainer(2,1);
  *zY << 0 << 1;
  vector<double> lambda(3);
  lambda[0] = 0;
  lambda[1] = 1;
  lambda[2] = 5;

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

  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;
  param.maxit          = 500;
  param.konv           = 0.0001;
  param.konvtime       = true;
  param.time           = false;
  param.test           = true;
  param.seconds        = 72000;
/*
  IsParameter paramgp;
  paramgp.lambdavalue    = 1.0;
  paramgp.lambdadeltaval = 1.0;
  paramgp.sigma          = 0.01;
  paramgp.maxit          = 500;
  paramgp.konv           = 0.0001;
  paramgp.time           = true;
  paramgp.test           = true;
  paramgp.seconds        = 60; */

  Test *testval     = new Test(2,2,10000,lambda, *zY,*zY,alphX, alphY);
  GIS *test         = new GIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
//  GISgp *testgp     = new GISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,paramgp);
  SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
//  SCGISgp *testscgp = new SCGISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,paramgp);

  for(int i=0;i<4;i++){
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->propAlphX(i,0)+test->propAlphX(i,1)+ test->propAlphX(i,2) +test->propAlphX(i,3),0.1);
 //       CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testgp->propAlphX(i,0)+testgp->propAlphX(i,1)+ testgp->propAlphX(i,2) + testgp->propAlphX(i,3),0.1);

        CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testsc->propAlphX(i,0)+testsc->propAlphX(i,1)+ testsc->propAlphX(i,2) +testsc->propAlphX(i,3),0.1);
 //       CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testscgp->propAlphX(i,0)+testscgp->propAlphX(i,1)+ testscgp->propAlphX(i,2) +testscgp->propAlphX(i,3),0.1);
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
  for(int i=0;i< 4; i++){
	  for(int j=0; j< 4; j++){
		  cout << i << " " << j << endl;
		  cout << "exakt: " <<  testval->getProp(i,j) << endl;
		  cout << "gis:  " << test->propAlphX(i,j) << endl;
		  cout << "scgis: " << testsc->propAlphX(i,j) << endl;
//		  cout << "gisgp: " << testgp->propAlphX(i,j) << endl;
//		  cout << "scgisgp: " << testscgp->propAlphX(i,j) << endl;
	  }
  }


  cout << " GIS: " << testval->KL(test);
  cout << " SCGIS: " << testval->KL(testsc);
//  cout << " GISgp: " << testval->KL(testgp);
//  cout << " SCGISgp: " << testval->KL(testscgp);
}

void itTest::NotBinarySix()
{
  cout << endl;
  cout << "NotBinary" << endl;
  srand(time(NULL));
  DContainer *zX = new DContainer(6,1);
  *zX << 0 << 1 << 2 << 3 << 4 << 5;
  DContainer *zY = new DContainer(6,1);
  *zY << 0 << 1 << 2 << 3 << 4 << 5;
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
  param.konvtime       = true;
  param.time           = false;
  param.test           = true;
  param.seconds        = 72000;

  Test *testval = new Test(1,1,10000,lambda, *zX,*zY,alphX, alphY);
  GIS *test     = new GIS(testval->getvalX(),testval->getvalY(),*zX,*zY,alphX,alphY,param);
  SCGIS *testsc = new SCGIS(testval->getvalX(),testval->getvalY(),*zX,*zY,alphX,alphY,param);

  for(int i=0;i<6;i++){
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1,test->propAlphX(i,0)+test->propAlphX(i,1)+test->propAlphX(i,2)+test->propAlphX(i,3)+test->propAlphX(i,4)+test->propAlphX(i,5),0.1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1,testsc->propAlphX(i,0)+testsc->propAlphX(i,1)+testsc->propAlphX(i,2)+testsc->propAlphX(i,3)+testsc->propAlphX(i,4)+testsc->propAlphX(i,5),0.1);
  }

  cout << " Wahrscheinlichkeiten: " << endl;
  for(int i=0;i< 6; i++){
	  for(int j=0; j< 6; j++){
		  cout << i << " " << j << endl;
		  cout << "exakt: " <<  testval->getProp(i,j) << endl;
		  cout << "gis:  " << test->propAlphX(i,j) << endl;
		  cout << "scgis: " << testsc->propAlphX(i,j) << endl;
	  }
  }



  for(int i=0; i< test->getsizeconv()-1; i++ )
  {
    CPPUNIT_ASSERT(test->getconv(i) >= test->getconv(i+1));
  }
  cout << " GIS: " << testval->KL(test);
  cout << " SCGIS: " << testval->KL(testsc);
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

  Test *testval     = new Test(4,4,100,lambda, *zX,*zY,alphX, alphY);
  GIS *test         = new GIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
 // GISgp *testgp     = new GISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
  SCGIS *testsc     = new SCGIS(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);
//  SCGISgp *testscgp = new SCGISgp(testval->getvalX(),testval->getvalY(),*zY,*zY,alphX,alphY,param);

  double sum1 = 0.0;
  double sum2 = 0.0;
  double sum3 = 0.0;
  double sum4 = 0.0;

      for(int k=0;k<16;k++){
    	  for(int i=0;i<81;i++){
   		    sum1+=test->propAlphX(k,i);
   		//    sum2+=testgp->propAlphX(k,i);
            sum3+=testsc->propAlphX(k,i);
   		//    sum4+=testscgp->propAlphX(k,i);
    	  }
          CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum1,0.1);
         // CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum2,0.1);
          CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum3,0.1);
        //  CPPUNIT_ASSERT_DOUBLES_EQUAL  (  1,sum4,0.1);
          sum1 = 0;
          sum2 = 0;
          sum3 = 0;
          sum4 = 0;
  }
  cout << " Wahrscheinlichkeiten: " << endl;
  for(int i=0;i< 6; i++){
     for(int j=0; j< 6; j++){
       cout << i << " " << j << endl;
       cout << "exakt: " <<  testval->getProp(i,j) << endl;
       cout << "gis:  " << test->propAlphX(i,j) << endl;
       cout << "scgis: " << testsc->propAlphX(i,j) << endl;
    	  }
      }

	      cout << " GIS: " << testval->KL(test);
	   //   cout << " GISgp: " << testval->KL(testgp);
	      cout << " SCGIS " << testval->KL(testsc);
	    //  cout << " SCGISgp " << testval->KL(testscgp);
}

*/
