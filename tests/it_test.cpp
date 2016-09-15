#include "it_test.h"

#include <entropy++/Container.h>
#include <entropy++/MI.h>
#include <entropy++/sparse/MI.h>
#include "../experiments/it/GIS.h"
#include "../experiments/it/SCGIS.h"
#include "../experiments/it/Comp.h"

#include <iostream>
#include <string>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( itTest );

void itTest::COMP(){

}
void itTest::OneXOneY()
{	cout << "OneXOneY" << endl;
	 srand(time(NULL));
	 DContainer *zX = new DContainer(2,1);
	 *zX << 0 << 1;
	  DContainer *zY = new DContainer(2,1);
	 *zY << 0 << 1;
	 vector<double> lambda(3);
	 lambda[0] = 0;
	 lambda[1] =1 ;
	 lambda[2] = 5;

	 Comp *zTest = new Comp(1,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,0);
	 /*
	 cout << "GIS "<< endl;
	 cout <<zTest->prop(0,0,0,0) << endl;
	 cout <<zTest->prop(0,0,1,0) << endl;
	 cout <<zTest->prop(0,0,0,1) << endl;
	 cout <<zTest->prop(0,0,1,1)<< endl;
	 cout << endl;
	 cout << endl;
	 cout <<zTest->getFeatureArraylambda(0,0,0,0) << endl;
	 cout <<zTest->getFeatureArraylambda(0,0,1,0) << endl;
	 cout <<zTest->getFeatureArraylambda(0,0,0,1) << endl;
	 cout <<zTest->getFeatureArraylambda(0,0,1,1) << endl;
	 cout << endl;
	*/
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,0,0,0)+zTest->prop(0,0,0,1,0),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,1,0,0)+zTest->prop(0,0,1,1,0),0.1);
	 for(int i=0; i< zTest->getsizeconv(0)-1; i++ ){
		if((zTest->getconv(i,0)>pow(10,-11))){
			if(!(zTest->getconv(i,0) >= zTest->getconv(i+1,0))){
				cout << zTest->getconv(i,0)  << " " << zTest->getconv(i+1,0) << endl;
			}
		CPPUNIT_ASSERT 	(zTest->getconv(i,0) >= zTest->getconv(i+1,0) );
		}
	 }
	 CPPUNIT_ASSERT (zTest->KL(0)<1.5);
}
void itTest::SCOneXOneY()
{cout << "SCOneXOneY" << endl;
	srand(time(NULL));
		 DContainer *zX = new DContainer(2,1);
		 *zX << 0 << 1;
		  DContainer *zY = new DContainer(2,1);
		 *zY << 0 << 1;
		 vector<double> lambda(3);
		 lambda[0] = 0;
		 lambda[1] =1 ;
		 lambda[2] = 5;

		 Comp *zTest = new Comp(1,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,1);
		 /*
		 cout << "SCGIS " << endl;
		 cout <<zTest->prop(0,0,0,0)-Test->prop(0,0,0,0) << endl;
		 cout <<zTest->prop(0,0,1,0)-Test->prop(0,0,1,0) << endl;
		 cout <<zTest->prop(0,0,0,1)-Test->prop(0,0,0,1) << endl;
		 cout <<zTest->prop(0,0,1,1)-Test->prop(0,0,1,1) << endl;
		 cout << endl;
		 cout << zTest->getFeatureArraylambda(0,0,0,0) <<endl;
		 cout << zTest->getFeatureArraylambda(0,0,1,0) <<endl;
		 cout << zTest->getFeatureArraylambda(0,0,0,1) <<endl;
		 cout << zTest->getFeatureArraylambda(0,0,1,1) <<endl;
		 cout << endl;
		*/
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,0,0,1)+zTest->prop(0,0,0,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,1,0,1)+zTest->prop(0,0,1,1,1),0.1);
;

		 for(int i=0; i< zTest->getsizeconv(1)-1; i++ ){
			if((zTest->getconv(i,1)>pow(10,-11))){
				if(!(zTest->getconv(i,1) >= zTest->getconv(i+1,1))){
					cout << zTest->getconv(i,1)  << " " << zTest->getconv(i+1,1) << endl;
				}
			CPPUNIT_ASSERT 	(zTest->getconv(i,1) >= zTest->getconv(i+1,1) );
			}
		 }
		 CPPUNIT_ASSERT (zTest->KL(1)<1.5);
}
void itTest::TwoXOneY()
{	cout << "TWoXOneY" << endl;
	 srand(time(NULL));
	 DContainer *zX = new DContainer(2,1);
	 *zX << 0 << 1;
	  DContainer *zY = new DContainer(2,1);
	 *zY << 0 << 1;
	 vector<double> lambda(3);
	 lambda[0] = 0;
	 lambda[1] =1 ;
	 lambda[2] = 5;

	 Comp *zTest = new Comp(2,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,0);

	 /*
	 cout << " GIS "<< endl;
	 cout <<zTest->prop(0,0,0,0)-Test->prop(0,0,0,0)  << endl;
	 cout <<zTest->prop(0,0,1,0)-Test->prop(0,0,1,0) << endl;
	 cout <<zTest->prop(0,0,0,1)-Test->prop(0,0,0,1) << endl;
	 cout <<zTest->prop(0,0,1,1)-Test->prop(0,0,1,1) << endl;
	 cout <<zTest->prop(1,0,0,0)-Test->prop(1,0,0,0) << endl;
	 cout <<zTest->prop(1,0,1,0)-Test->prop(1,0,1,0) << endl;
	 cout <<zTest->prop(1,0,0,1)-Test->prop(1,0,0,1) << endl;
	 cout <<zTest->prop(1,0,1,1)-Test->prop(1,0,1,1) << endl;
	 cout << endl;
	*/
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,0,0,0)+zTest->prop(0,0,0,1,0),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,1,0,0)+zTest->prop(0,0,1,1,0),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,0,0,0,0)+zTest->prop(1,0,0,1,0),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,0,1,0,0)+zTest->prop(1,0,1,1,0),0.1);

	 for(int i=0; i< zTest->getsizeconv(0)-1; i++ ){
		if((zTest->getconv(i,0)>pow(10,-11))){
			if(!(zTest->getconv(i,0) >= zTest->getconv(i+1,0))){
				cout << zTest->getconv(i,0)  << " " << zTest->getconv(i+1,0) << endl;
			}
		CPPUNIT_ASSERT 	(zTest->getconv(i,0) >= zTest->getconv(i+1,0) );
		}
	 }
	 CPPUNIT_ASSERT (zTest->KL(0)<1.5);
}
void itTest::SCTwoXOneY()
{	 cout << "SCTWoXOneY" << endl;
	srand(time(NULL));
	DContainer *zX = new DContainer(2,1);
		*zX << 0 << 1;
	DContainer *zY = new DContainer(2,1);
		*zY << 0 << 1;
	vector<double> lambda(3);
		lambda[0] = 0;
		lambda[1] =1 ;
		lambda[2] = 5;

	Comp *zTest = new Comp(2,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,1);

	 /*
	 cout << "SCGIS " <<  endl;
	 cout <<zTest->prop(0,0,0,0)-Test->prop(0,0,0,0)  << endl;
	 cout <<zTest->prop(0,0,1,0)-Test->prop(0,0,1,0) << endl;
	 cout <<zTest->prop(0,0,0,1)-Test->prop(0,0,0,1) << endl;
	 cout <<zTest->prop(0,0,1,1)-Test->prop(0,0,1,1) << endl;
	 cout <<zTest->prop(1,0,0,0)-Test->prop(1,0,0,0) << endl;
	 cout <<zTest->prop(1,0,1,0)-Test->prop(1,0,1,0) << endl;
	 cout <<zTest->prop(1,0,0,1)-Test->prop(1,0,0,1) << endl;
	 cout <<zTest->prop(1,0,1,1)-Test->prop(1,0,1,1) << endl;
	 cout << endl;
	*/
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,0,0,1)+zTest->prop(0,0,0,1,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,1,0,1)+zTest->prop(0,0,1,1,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,0,0,0,1)+zTest->prop(1,0,0,1,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,0,1,0,1)+zTest->prop(1,0,1,1,1),0.1);

	 for(int i=0; i< zTest->getsizeconv(1)-1; i++ ){
		if((zTest->getconv(i,1)>pow(10,-11))){
			if(!(zTest->getconv(i,1) >= zTest->getconv(i+1,1))){
				cout << zTest->getconv(i,1)  << " " << zTest->getconv(i+1,1) << endl;
			}
		CPPUNIT_ASSERT 	(zTest->getconv(i,1) >= zTest->getconv(i+1,1) );
		}
	 }
	 CPPUNIT_ASSERT (zTest->KL(1)<1.5);
}
void itTest::TwoXTwoY(){
	 cout << "TwoXTwoY" << endl;
	 srand(time(NULL));
	 DContainer *zX = new DContainer(2,1);
	 *zX << 0 << 1;
	  DContainer *zY = new DContainer(2,1);
	 *zY << 0 << 1;
	 vector<double> lambda(3);
	 lambda[0] = 0;
	 lambda[1] =1 ;
	 lambda[2] = 5;

	 Comp *zTest = new Comp(2,100,2,lambda,*zX,*zY,500,0.0001,false,true,0,0);


		 /*
			 cout << "GIS " << endl;
			 cout <<zTest->prop(0,0,0,0)-Test->prop(0,0,0,0)  << endl;
			 cout <<zTest->prop(0,0,1,0)-Test->prop(0,0,1,0) << endl;
			 cout <<zTest->prop(0,0,0,1)-Test->prop(0,0,0,1) << endl;
			 cout <<zTest->prop(0,0,1,1)-Test->prop(0,0,1,1) << endl;
			 cout << endl;
			 cout <<zTest->prop(1,0,0,0)-Test->prop(1,0,0,0) << endl;
			 cout <<zTest->prop(1,0,1,0)-Test->prop(1,0,1,0) << endl;
			 cout <<zTest->prop(1,0,0,1)-Test->prop(1,0,0,1) << endl;
			 cout <<zTest->prop(1,0,1,1)-Test->prop(1,0,1,1) << endl;
			 cout << endl;
			 cout <<zTest->prop(0,1,0,0)-Test->prop(0,1,0,0) << endl;
			 cout <<zTest->prop(0,1,1,0)-Test->prop(0,1,1,0) << endl;
			 cout <<zTest->prop(0,1,0,1)-Test->prop(0,1,0,1) << endl;
			 cout <<zTest->prop(0,1,1,1)-Test->prop(0,1,1,1) << endl;
			 cout << endl;
			 cout <<zTest->prop(1,1,0,0)-Test->prop(1,1,0,0) << endl;
			 cout <<zTest->prop(1,1,1,0)-Test->prop(1,1,1,0) << endl;
			 cout <<zTest->prop(1,1,0,1)-Test->prop(1,1,0,1) << endl;
			 cout <<zTest->prop(1,1,1,1)-Test->prop(1,1,1,1) << endl;
			 cout << endl;
			*/
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,0,0,0)+zTest->prop(0,0,0,1,0),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,1,0,0)+zTest->prop(0,0,1,1,0),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,0,0,0,0)+zTest->prop(1,0,0,1,0),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,0,1,0,0)+zTest->prop(1,0,1,1,0),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,1,0,0,0)+zTest->prop(0,1,0,1,0),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,1,1,0,0)+zTest->prop(0,1,1,1,0),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,1,0,0,0)+zTest->prop(1,1,0,1,0),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,1,1,0,0)+zTest->prop(1,1,1,1,0),0.1);


			 for(int i=0; i< zTest->getsizeconv(0)-1; i++ ){
				if((zTest->getconv(i,0)>pow(10,-11))){
					if(!(zTest->getconv(i,0) >= zTest->getconv(i+1,0))){
						cout << zTest->getconv(i,0)  << " " << zTest->getconv(i+1,0) << endl;
					}
				CPPUNIT_ASSERT 	(zTest->getconv(i,0) >= zTest->getconv(i+1,0) );
				}
			 }
			 CPPUNIT_ASSERT (zTest->KL(0)<1.5);
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
	 lambda[1] =1 ;
	 lambda[2] = 5;

	 Comp *zTest = new Comp(2,100,2,lambda,*zX,*zY,500,0.0001,false,true,0,1);
	 /*
		 cout << "SCGIS " << endl;
		 cout <<zTest->prop(0,0,0,0)-Test->prop(0,0,0,0)  << endl;
		 cout <<zTest->prop(0,0,1,0)-Test->prop(0,0,1,0) << endl;
		 cout <<zTest->prop(0,0,0,1)-Test->prop(0,0,0,1) << endl;
		 cout <<zTest->prop(0,0,1,1)-Test->prop(0,0,1,1) << endl;
		 cout << endl;
		 cout <<zTest->prop(1,0,0,0)-Test->prop(1,0,0,0) << endl;
		 cout <<zTest->prop(1,0,1,0)-Test->prop(1,0,1,0) << endl;
		 cout <<zTest->prop(1,0,0,1)-Test->prop(1,0,0,1) << endl;
		 cout <<zTest->prop(1,0,1,1)-Test->prop(1,0,1,1) << endl;
		 cout << endl;
		 cout <<zTest->prop(0,1,0,0)-Test->prop(0,1,0,0) << endl;
		 cout <<zTest->prop(0,1,1,0)-Test->prop(0,1,1,0) << endl;
		 cout <<zTest->prop(0,1,0,1)-Test->prop(0,1,0,1) << endl;
		 cout <<zTest->prop(0,1,1,1)-Test->prop(0,1,1,1) << endl;
		 cout << endl;
		 cout <<zTest->prop(1,1,0,0)-Test->prop(1,1,0,0) << endl;
		 cout <<zTest->prop(1,1,1,0)-Test->prop(1,1,1,0) << endl;
		 cout <<zTest->prop(1,1,0,1)-Test->prop(1,1,0,1) << endl;
		 cout <<zTest->prop(1,1,1,1)-Test->prop(1,1,1,1) << endl;
		 cout << endl;
		 */
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,0,0,1)+zTest->prop(0,0,0,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,1,0,1)+zTest->prop(0,0,1,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,0,0,0,1)+zTest->prop(1,0,0,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,0,1,0,1)+zTest->prop(1,0,1,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,1,0,0,1)+zTest->prop(0,1,0,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,1,1,0,1)+zTest->prop(0,1,1,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,1,0,0,1)+zTest->prop(1,1,0,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(1,1,1,0,1)+zTest->prop(1,1,1,1,1),0.1);

		 for(int i=0; i< zTest->getsizeconv(1)-1; i++ ){
			if((zTest->getconv(i,1)> pow(10,-11))){
				if(!(zTest->getconv(i,1) >= zTest->getconv(i+1,1))){
					cout << zTest->getconv(i,1)  << " " << zTest->getconv(i+1,1) << endl;
				}
			CPPUNIT_ASSERT 	(zTest->getconv(i,1) >= zTest->getconv(i+1,1) );
			}
		 }
		 CPPUNIT_ASSERT (zTest->KL(1)<1.5);

}
void itTest::NotBinary(){
	cout << "NotBinary" << endl;
	srand(time(NULL));
	 DContainer *zX = new DContainer(4,1);
	 *zX << 0 << 1 << 2 << 3;
	  DContainer *zY = new DContainer(2,1);
	 *zY << 0 << 1;
	 vector<double> lambda(3);
	 lambda[0] = 0;
	 lambda[1] =1 ;
	 lambda[2] = 5;

	 Comp *zTest = new Comp(1,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,0);

		 //cout << endl;
		 //cout <<zTest->gis(0,0,0,0)-Test->gis(0,0,0,0) << endl;
		 //cout <<zTest->gis(0,0,0,1)-Test->gis(0,0,0,1) << endl;
		 //cout <<zTest->gis(0,0,1,0)-Test->gis(0,0,1,0) << endl;
		 //cout <<zTest->gis(0,0,1,1)-Test->gis(0,0,1,1) << endl;
		 //cout <<zTest->gis(0,0,2,0)-Test->gis(0,0,0,0) << endl;
		 //cout <<zTest->gis(0,0,2,1)-Test->gis(0,0,0,1) << endl;
		 //cout <<zTest->gis(0,0,3,0)-Test->gis(0,0,1,0) << endl;
		 //cout <<zTest->gis(0,0,3,1)-Test->gis(0,0,1,1) << endl;
		 //cout << endl;

		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,0,0,0)+zTest->prop(0,0,0,1,0),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,1,0,0)+zTest->prop(0,0,1,1,0),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,2,0,0)+zTest->prop(0,0,2,1,0),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->prop(0,0,3,0,0)+zTest->prop(0,0,3,1,0),0.1);


		 for(int i=0; i< zTest->getsizeconv(0)-1; i++ ){
			 CPPUNIT_ASSERT 	(zTest->getconv(i,0)>= zTest->getconv(i+1,0) );
		 }
		 CPPUNIT_ASSERT (zTest->KL(0)<1.5);
	}
