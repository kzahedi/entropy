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

	 Comp *Test = new Comp(1,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,0);
	 Comp *Testgp= new Comp(1,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,2);

	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Test->prop(0,0,(*zX)(0,0),(*zY)(0,0))+Test->prop(0,0,(*zX)(0,0),(*zY)(1,0)),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Test->prop(0,0,(*zX)(1,0),(*zY)(0,0))+Test->prop(0,0,(*zX)(1,0),(*zY)(1,0)),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testgp->prop(0,0,(*zX)(0,0),(*zY)(0,0))+Testgp->prop(0,0,(*zX)(0,0),(*zY)(1,0)),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testgp->prop(0,0,(*zX)(1,0),(*zY)(0,0))+Testgp->prop(0,0,(*zX)(1,0),(*zY)(1,0)),0.1);
	 for(int i=0; i< Test->getsizeconv()-1; i++ ){
		if((Test->getconv(i)>pow(10,-11))){
			if(!(Test->getconv(i) >= Test->getconv(i+1))){
				cout << Test->getconv(i)  << " " << Test->getconv(i+1) << endl;
			}
		CPPUNIT_ASSERT 	(Test->getconv(i) >= Test->getconv(i+1) );
		}
	 }
	 CPPUNIT_ASSERT (Test->KL1()<1.5);
	 CPPUNIT_ASSERT (Testgp->KL1()<2);
}
void itTest::SCOneXOneY()
{ cout << "SCOneXOneY" << endl;
	srand(time(NULL));
		 DContainer *zX = new DContainer(2,1);
		 *zX << 0 << 1;
		  DContainer *zY = new DContainer(2,1);
		 *zY << 0 << 1;
		 vector<double> lambda(3);
		 lambda[0] = 0;
		 lambda[1] =1 ;
		 lambda[2] = 5;

		 Comp *Test = new Comp(1,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,1);
		 Comp *Testgp = new Comp(1,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,3);

		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Test->prop(0,0,(*zX)(0,0),(*zY)(0,0))+Test->prop(0,0,(*zX)(0,0),(*zY)(1,0)),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Test->prop(0,0,(*zX)(1,0),(*zY)(0,0))+Test->prop(0,0,(*zX)(1,0),(*zY)(1,0)),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testgp->prop(0,0,(*zX)(0,0),(*zY)(0,0))+Testgp->prop(0,0,(*zX)(0,0),(*zY)(1,0)),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testgp->prop(0,0,(*zX)(1,0),(*zY)(0,0))+Testgp->prop(0,0,(*zX)(1,0),(*zY)(1,0)),0.1);

		 for(int i=0; i< Test->getsizeconv()-1; i++ ){
			if((Test->getconv(i)>pow(10,-11))){
				if(!(Test->getconv(i) >= Test->getconv(i+1))){
					cout << Test->getconv(i)  << " " << Test->getconv(i+1) << endl;
				}
			CPPUNIT_ASSERT 	(Test->getconv(i) >= Test->getconv(i+1) );
			}
		 }
		 CPPUNIT_ASSERT (Test->KL1()<1.5);
		 CPPUNIT_ASSERT (Testgp->KL1()<2);
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

	 Comp *Test = new Comp(2,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,0);
	 Comp *Testgp = new Comp(2,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,2);

	 for(int i=0;i<2;i++){
			 for(int k=0;k<2;k++){
					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Test->prop(i,0,(*zX)(k,0),(*zY)(0,0))+Test->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testgp->prop(i,0,(*zX)(k,0),(*zY)(0,0))+Testgp->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);

			 }
	 }

	 for(int i=0; i< Test->getsizeconv()-1; i++ ){
		if((Test->getconv(i)>pow(10,-11))){
			if(!(Test->getconv(i) >= Test->getconv(i+1))){
				cout << Test->getconv(i)  << " " << Test->getconv(i+1) << endl;
			}
		CPPUNIT_ASSERT 	(Test->getconv(i) >= Test->getconv(i+1) );
		}
	 }
	 CPPUNIT_ASSERT (Test->KL1()<1.5);
	 CPPUNIT_ASSERT (Testgp->KL1()<2);
}
void itTest::SCTwoXOneY()
{ 	 cout << "SCTWoXOneY" << endl;
	srand(time(NULL));
	DContainer *zX = new DContainer(2,1);
		*zX << 0 << 1;
	DContainer *zY = new DContainer(2,1);
		*zY << 0 << 1;
	vector<double> lambda(3);
		lambda[0] = 0;
		lambda[1] =1 ;
		lambda[2] = 5;

	Comp *Test = new Comp(2,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,1);
	Comp *Testgp = new Comp(2,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,3);

	 for(int i=0;i<2;i++){
			 for(int k=0;k<2;k++){
					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Test->prop(i,0,(*zX)(k,0),(*zY)(0,0))+Test->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testgp->prop(i,0,(*zX)(k,0),(*zY)(0,0))+Testgp->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);

			 }
	 }

	 for(int i=0; i< Test->getsizeconv()-1; i++ ){
		if((Test->getconv(i)>pow(10,-11))){
			if(!(Test->getconv(i) >= Test->getconv(i+1))){
				cout << Test->getconv(i)  << " " << Test->getconv(i+1) << endl;
			}
		CPPUNIT_ASSERT 	(Test->getconv(i) >= Test->getconv(i+1) );
		}
	 }
	 CPPUNIT_ASSERT (Test->KL1()<2);
	 CPPUNIT_ASSERT (Testgp->KL1()<2);

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

	 Comp *Test = new Comp(2,100,2,lambda,*zX,*zY,500,0.0001,false,true,0,0);
	 Comp *Testgp = new Comp(2,100,2,lambda,*zX,*zY,500,0.0001,false,true,0,2);

	 	 	 for(int i=0;i<2;i++){
	 	 		 for(int j=0;j<2;j++){
	 	 			 for(int k=0;k<2;k++){
	 					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Test->prop(i,0,(*zX)(k,0),(*zY)(0,0))+Test->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
	 					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testgp->prop(i,0,(*zX)(k,0),(*zY)(0,0))+Testgp->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
	 	 			 }
	 	 		 }
	 	 	 }

			 for(int i=0; i< Test->getsizeconv()-1; i++ ){
				if((Test->getconv(i)>pow(10,-11))){
					if(!(Test->getconv(i) >= Test->getconv(i+1))){
						cout << Test->getconv(i)  << " " << Test->getconv(i+1) << endl;
					}
				CPPUNIT_ASSERT 	(Test->getconv(i) >= Test->getconv(i+1) );
				}
			 }
			 CPPUNIT_ASSERT (Test->KL1()<2);
			 CPPUNIT_ASSERT (Testgp->KL1()<2);

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

	 Comp *Test = new Comp(2,100,2,lambda,*zX,*zY,500,0.0001,false,true,0,1);
	 Comp *Testgp = new Comp(2,100,2,lambda,*zX,*zY,500,0.0001,false,true,0,3);


	 for(int i=0;i<2;i++){
		 for(int j=0;j<2;j++){
			 for(int k=0;k<2;k++){
					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Test->prop(i,0,(*zX)(k,0),(*zY)(0,0))+Test->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);
					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testgp->prop(i,0,(*zX)(k,0),(*zY)(0,0))+Testgp->prop(i,0,(*zX)(k,0),(*zY)(1,0)),0.1);

			 }
		 }
	 }

		 for(int i=0; i< Test->getsizeconv()-1; i++ ){
			if((Test->getconv(i)> pow(10,-11))){
				if(!(Test->getconv(i) >= Test->getconv(i+1))){
					cout << Test->getconv(i)  << " " << Test->getconv(i+1) << endl;
				}
			CPPUNIT_ASSERT 	(Test->getconv(i) >= Test->getconv(i+1) );
			}
		 }
		 CPPUNIT_ASSERT (Test->KL1()<2);
		 CPPUNIT_ASSERT (Testgp->KL1()<2);

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

	 Comp *Test = new Comp(1,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,0);
	 Comp *Testgp = new Comp(1,100,1,lambda,*zX,*zY,500,0.0001,false,true,0,2);

	 for(int i=0;i<4;i++){
			CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Test->prop(0,0,(*zX)(i,0),(*zY)(0,0))+Test->prop(0,0,(*zX)(i,0),(*zY)(1,0)),0.1);
			CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testgp->prop(0,0,(*zX)(i,0),(*zY)(0,0))+Testgp->prop(0,0,(*zX)(i,0),(*zY)(1,0)),0.1);
	 }


		 for(int i=0; i< Test->getsizeconv()-1; i++ ){
			 CPPUNIT_ASSERT 	(Test->getconv(i)>= Test->getconv(i+1) );
		 }
		 CPPUNIT_ASSERT (Test->KL1()<1.5);
		 CPPUNIT_ASSERT (Testgp->KL1()<1.5);

	}
void itTest::FourXFourY(){

	cout << "FourXFourY" << endl;
	 srand(time(NULL));
	 DContainer *zX = new DContainer(2,1);
	 *zX << 0 << 1;
	  DContainer *zY = new DContainer(3,1);
	 *zY << 0 << 1 << 2;
	 vector<double> lambda(3);
	 lambda[0] = 0;
	 lambda[1] =1 ;
	 lambda[2] = 5;

	 Comp *Test = new Comp(4,100,4,lambda,*zX,*zY,200,0.0001,false,true,0,0);
	 Comp *Testsc = new Comp(4,100,4,lambda,*zX,*zY,200,0.0001,false,true,0,1);
	 Comp *Testgp = new Comp(4,100,4,lambda,*zX,*zY,200,0.0001,false,true,0,2);
	 Comp *Testscgp = new Comp(4,100,4,lambda,*zX,*zY,200,0.0001,false,true,0,3);

	 for(int i=0;i<4;i++){
		 for(int j=0;j<4;j++){
			 for(int k=0;k<2;k++){
					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Test->prop(i,j,k,0)+Test->prop(i,j,k,1)+Test->prop(i,j,k,2),0.1);
					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testsc->prop(i,j,k,0)+Testsc->prop(i,j,k,1)+Testsc->prop(i,j,k,2),0.1);
					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testgp->prop(i,j,k,0)+Testgp->prop(i,j,k,1)+Testgp->prop(i,j,k,2),0.1);
					CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,Testscgp->prop(i,j,k,0)+Testscgp->prop(i,j,k,1)+Testscgp->prop(i,j,k,2),0.1);

			 }
		 }
	 }
		 CPPUNIT_ASSERT (Test->KL1()<2);
		 CPPUNIT_ASSERT (Testsc->KL1()<2);
		 CPPUNIT_ASSERT (Testgp->KL1()<2);
		 CPPUNIT_ASSERT (Testscgp->KL1()<2);
}
