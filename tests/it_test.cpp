#include "it_test.h"

#include <entropy++/Container.h>
#include <entropy++/MI.h>
#include <entropy++/sparse/MI.h>
#include "../experiments/it/GIS.h"

#include <iostream>
#include <string>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( itTest );


void itTest:: OneXOneY()
{
	 srand(time(NULL));
	 int n=100;
	 DContainer *eX = new DContainer(n,1);
	 for(int i=0;i< n;i++ ){
		 for(int j=0;j<1;j++){
			 *eX << rand() % 2;
		 }
	 }
	 DContainer *zX = new DContainer(2,1);
	 *zX << 0 << 1;
	  DContainer *zY = new DContainer(2,1);
	 *zY << 0 << 1;

	 GIS *Test = new GIS(1,*eX);
	 Test->setFeatureArraylambda(0,0,1,0,0);
	 Test->setFeatureArraylambda(0,0,1,1,0);
	 Test->setFeatureArraylambda(0,0,0,0,1);
	 Test->setFeatureArraylambda(0,0,0,1,1);

	 double** prop;
	 prop=new double*[n];
	 for(int i=0;i<n;i++){
		 prop[i]=new double[2];
		 for(int j=0;j<2;j++){
			 prop[i][j]=0;
		 }
	 }

	 for(int i=0;i<n;i++){
		for(int propi=0;propi<2;propi++){
				prop[i][propi]=Test->gis(0,0,(*eX)(i,0),propi);
			 }
		 }
	 DContainer *esY=new DContainer(n,1);

	 for(int i=0;i<n;i++){
		 double z=(double)rand()/RAND_MAX;
		 double s=0;
		 int ind=0;
		 for(int j=0;j<2 && s<z;j++){
				s+=prop[i][j];
				ind=j;
			 }
		 (*esY) << ind;
		 }

	 GIS *zTest = new GIS(*eX,*esY,*zX,*zY,1,5000,0.01);

	 cout << endl;
	 cout <<zTest->gis(0,0,0,0)-Test->gis(0,0,0,0)  << endl;
	 cout <<zTest->gis(0,0,1,0)-Test->gis(0,0,1,0) << endl;
	 cout <<zTest->gis(0,0,0,1)-Test->gis(0,0,0,1)  << endl;
	 cout <<zTest->gis(0,0,1,1)-Test->gis(0,0,1,1) << endl;
	 cout << endl;
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,0,0)+zTest->gis(0,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,1,0)+zTest->gis(0,0,1,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,0,0) ,zTest->gis(0,0,0,0),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,0,1) ,zTest->gis(0,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,1,0) ,zTest->gis(0,0,1,0),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,1,1) ,zTest->gis(0,0,1,1),0.1);

}

void itTest:: TwoXOneY()
{
	 srand(time(NULL));
	 int n=100;
	 DContainer *eX = new DContainer(n,2);
	 for(int i=0;i< n;i++ ){
		 for(int j=0;j<2;j++){
			 *eX << rand() % 2;
		 }
	 }

	 DContainer *zX = new DContainer(2,1);
	 *zX << 0 << 1;
	  DContainer *zY = new DContainer(2,1);
	 *zY << 0 << 1;

	 GIS *Test = new GIS(1,*eX);

	 Test->setFeatureArraylambda(0,0,1,0,4);
	 Test->setFeatureArraylambda(0,0,1,1,0);
	 Test->setFeatureArraylambda(0,0,0,0,1);
	 Test->setFeatureArraylambda(0,0,0,1,5);
	 Test->setFeatureArraylambda(1,0,1,0,2);
	 Test->setFeatureArraylambda(1,0,1,1,3);
	 Test->setFeatureArraylambda(1,0,0,0,0);
	 Test->setFeatureArraylambda(1,0,0,1,1);

	 double** prop;
	 prop=new double*[n];
	 for(int i=0;i<n;i++){
		 prop[i]=new double[2];
		 for(int j=0;j<2;j++){
			 prop[i][j]=0;
		 }
	 }
	 for(int i=0;i<n;i++){
		for(int propi=0;propi<2;propi++){
				prop[i][propi]=Test->gis(0,0,(*eX)(i,0),propi)*Test->gis(1,0,(*eX)(i,1),propi);
		}
	 }
	 DContainer *esY=new DContainer(n,1);

	 for(int i=0;i<n;i++){
		 double z=(double)rand()/RAND_MAX;
		 double s=0;
		 int ind=0;
		 for(int j=0;j<2 && s<z;j++){
				s+=prop[i][j];
				ind=j;
		 }
		 (*esY) << ind;

	 }
	 GIS *zTest = new GIS(*eX,*esY,*zX,*zY,1,5000,0.01);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,0,0)+zTest->gis(0,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,1,0)+zTest->gis(0,0,1,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(1,0,0,0)+zTest->gis(1,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(1,0,1,0)+zTest->gis(1,0,1,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,0,0) ,zTest->gis(0,0,0,0),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,0,1) ,zTest->gis(0,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,1,0) ,zTest->gis(0,0,1,0),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,1,1) ,zTest->gis(0,0,1,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(1,0,0,0) ,zTest->gis(1,0,0,0),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(1,0,0,1) ,zTest->gis(1,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(1,0,1,0) ,zTest->gis(1,0,1,0),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(1,0,1,1) ,zTest->gis(1,0,1,1),0.1);

	 cout << endl;
	 cout <<zTest->gis(0,0,0,0)-Test->gis(0,0,0,0)  << endl;
	 cout <<zTest->gis(0,0,1,0)-Test->gis(0,0,1,0) << endl;
	 cout <<zTest->gis(0,0,0,1)-Test->gis(0,0,0,1) << endl;
	 cout <<zTest->gis(0,0,1,1)-Test->gis(0,0,1,1) << endl;
	 cout <<zTest->gis(1,0,0,0)-Test->gis(1,0,0,0) << endl;
	 cout <<zTest->gis(1,0,1,0)-Test->gis(1,0,1,0) << endl;
	 cout <<zTest->gis(1,0,0,1)-Test->gis(1,0,0,1) << endl;
	 cout <<zTest->gis(1,0,1,1)-Test->gis(1,0,1,1) << endl;
	 cout << endl;




}

