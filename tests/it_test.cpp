#include "it_test.h"

#include <entropy++/Container.h>
#include <entropy++/MI.h>
#include <entropy++/sparse/MI.h>
#include "../experiments/it/GIS.h"
#include "../experiments/it/SCGIS.h"

#include <iostream>
#include <string>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( itTest );


void itTest::OneXOneY()
{
	 srand(time(NULL));
	 int n=1000;
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

	 GIS *Test = new GIS(1,*eX,*zX,*zY);
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
	 GIS *zTest = new GIS(*eX,*esY,*zX,*zY,1,500,0.01,true);

	 //cout << endl;
	 //cout <<zTest->gis(0,0,0,0)-Test->gis(0,0,0,0) << endl;
	 //cout <<zTest->gis(0,0,1,0)-Test->gis(0,0,1,0) << endl;
	 //cout <<zTest->gis(0,0,0,1)-Test->gis(0,0,0,1) << endl;
	 //cout <<zTest->gis(0,0,1,1)-Test->gis(0,0,1,1) << endl;
	 //cout << endl;

	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,0,0)+zTest->gis(0,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,1,0)+zTest->gis(0,0,1,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,0,0) ,zTest->gis(0,0,0,0),0.2);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,0,1) ,zTest->gis(0,0,0,1),0.2);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,1,0) ,zTest->gis(0,0,1,0),0.2);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,1,1) ,zTest->gis(0,0,1,1),0.2);

	 for(int i=0; i< zTest->getsizeconv()-1; i++ ){
		 CPPUNIT_ASSERT 	(zTest->getconv(i)>= zTest->getconv(i+1) );
	 }
	 Test->~GIS();
	 zTest->~GIS();
}
void itTest::SCOneXOneY()
{
	srand(time(NULL));
		 int n=1000;
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

		 GIS *Test = new GIS(1,*eX,*zX,*zY);
		 Test->setFeatureArraylambda(0,0,1,0,2);
		 Test->setFeatureArraylambda(0,0,1,1,0);
		 Test->setFeatureArraylambda(0,0,0,0,1);
		 Test->setFeatureArraylambda(0,0,0,1,3);

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
		 SCGIS *zTest = new SCGIS(*eX,*esY,*zX,*zY,1,5000,0.01,true);

		 cout << endl;
		 cout <<zTest->scgis(0,0,0,0)-Test->gis(0,0,0,0) << endl;
		 cout <<zTest->scgis(0,0,1,0)-Test->gis(0,0,1,0) << endl;
		 cout <<zTest->scgis(0,0,0,1)-Test->gis(0,0,0,1) << endl;
		 cout <<zTest->scgis(0,0,1,1)-Test->gis(0,0,1,1) << endl;
		 cout << endl;
		 //cout << zTest->getFeatureArraylambda(0,0,0,0) <<endl;
		 //cout << zTest->getFeatureArraylambda(0,0,1,0) <<endl;
		 //cout << zTest->getFeatureArraylambda(0,0,0,1) <<endl;
		 //cout << zTest->getFeatureArraylambda(0,0,1,1) <<endl;
		 //cout << endl;

		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(0,0,0,0)+zTest->scgis(0,0,0,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(0,0,1,0)+zTest->scgis(0,0,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,0,0) ,zTest->scgis(0,0,0,0),0.2);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,0,1) ,zTest->scgis(0,0,0,1),0.2);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,1,0) ,zTest->scgis(0,0,1,0),0.2);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,1,1) ,zTest->scgis(0,0,1,1),0.2);

		 for(int i=0; i< zTest->getsizeconv()-1; i++ ){
			 CPPUNIT_ASSERT 	(zTest->getconv(i)>= zTest->getconv(i+1) );
		 }
		 Test->~GIS();
		 zTest->~SCGIS();
}
void itTest::TwoXOneY()
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

	 GIS *Test = new GIS(1,*eX,*zX,*zY);

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
	 vector<vector<double> > val(2,vector<double>(1));
	 val[0][0]=0;
	 val[1][0]=1;
	 for(int i=0;i<n;i++){
		for(int propi=0;propi<2;propi++){
				prop[i][propi]=Test->gis(i,val,propi);
		}
		CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,prop[i][0]+prop[i][1],0.1);
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
	 GIS *zTest = new GIS(*eX,*esY,*zX,*zY,1,500,0.01,true);
	 //cout << endl;
	 //cout <<zTest->gis(0,0,0,0)-Test->gis(0,0,0,0)  << endl;
	 //cout <<zTest->gis(0,0,1,0)-Test->gis(0,0,1,0) << endl;
	 //cout <<zTest->gis(0,0,0,1)-Test->gis(0,0,0,1) << endl;
	 //cout <<zTest->gis(0,0,1,1)-Test->gis(0,0,1,1) << endl;
	 //cout <<zTest->gis(1,0,0,0)-Test->gis(1,0,0,0) << endl;
	 //cout <<zTest->gis(1,0,1,0)-Test->gis(1,0,1,0) << endl;
	 //cout <<zTest->gis(1,0,0,1)-Test->gis(1,0,0,1) << endl;
	 //cout <<zTest->gis(1,0,1,1)-Test->gis(1,0,1,1) << endl;
	 //cout << endl;

	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,0,0)+zTest->gis(0,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,1,0)+zTest->gis(0,0,1,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(1,0,0,0)+zTest->gis(1,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(1,0,1,0)+zTest->gis(1,0,1,1),0.1);

	 for(int i=0; i< zTest->getsizeconv()-1; i++ ){
		 CPPUNIT_ASSERT 	(zTest->getconv(i)>= zTest->getconv(i+1) );
	 }
	 Test->~GIS();
	 zTest->~GIS();
}
void itTest::SCTwoXOneY()
{	 srand(time(NULL));
	 int n=1000;
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

	 GIS *Test = new GIS(1,*eX,*zX,*zY);

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
	 vector<vector<double> > val(2,vector<double>(1));
	 val[0][0]=0;
	 val[1][0]=1;
	 for(int i=0;i<n;i++){
		for(int propi=0;propi<2;propi++){
				prop[i][propi]=Test->gis(i,val,propi);
		}
		CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,prop[i][0]+prop[i][1],0.1);
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
	 SCGIS *zTest = new SCGIS(*eX,*esY,*zX,*zY,1,5000,0.01,true);
	 cout << endl;
	 cout <<zTest->scgis(0,0,0,0)-Test->gis(0,0,0,0)  << endl;
	 cout <<zTest->scgis(0,0,1,0)-Test->gis(0,0,1,0) << endl;
	 cout <<zTest->scgis(0,0,0,1)-Test->gis(0,0,0,1) << endl;
	 cout <<zTest->scgis(0,0,1,1)-Test->gis(0,0,1,1) << endl;
	 cout <<zTest->scgis(1,0,0,0)-Test->gis(1,0,0,0) << endl;
	 cout <<zTest->scgis(1,0,1,0)-Test->gis(1,0,1,0) << endl;
	 cout <<zTest->scgis(1,0,0,1)-Test->gis(1,0,0,1) << endl;
	 cout <<zTest->scgis(1,0,1,1)-Test->gis(1,0,1,1) << endl;
	 cout << endl;

	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(0,0,0,0)+zTest->scgis(0,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(0,0,1,0)+zTest->scgis(0,0,1,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(1,0,0,0)+zTest->scgis(1,0,0,1),0.1);
	 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(1,0,1,0)+zTest->scgis(1,0,1,1),0.1);

	 for(int i=0; i< zTest->getsizeconv()-1; i++ ){
		// cout << "conv "<< zTest->getconv(i) << " " << zTest->getconv(i+1) << endl;
		 CPPUNIT_ASSERT 	(zTest->getconv(i)>= zTest->getconv(i+1) );
	 }
	 Test->~GIS();
	 zTest->~SCGIS();

}
void itTest::TwoXTwoY(){
	 srand(time(NULL));
	 int n=1000;
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

	 GIS *Test = new GIS(2,*eX,*zX,*zY);

		 Test->setFeatureArraylambda(0,0,1,0,4);
		 Test->setFeatureArraylambda(0,0,1,1,0);
		 Test->setFeatureArraylambda(0,0,0,0,1);
		 Test->setFeatureArraylambda(0,0,0,1,5);

		 Test->setFeatureArraylambda(1,0,1,0,2);
		 Test->setFeatureArraylambda(1,0,1,1,3);
		 Test->setFeatureArraylambda(1,0,0,0,0);
		 Test->setFeatureArraylambda(1,0,0,1,1);

		 Test->setFeatureArraylambda(0,1,1,0,0);
		 Test->setFeatureArraylambda(0,1,1,1,0);
		 Test->setFeatureArraylambda(0,1,0,0,3);
		 Test->setFeatureArraylambda(0,1,0,1,5);

		 Test->setFeatureArraylambda(1,1,1,0,1);
		 Test->setFeatureArraylambda(1,1,1,1,3);
		 Test->setFeatureArraylambda(1,1,0,0,2);
		 Test->setFeatureArraylambda(1,1,0,1,1);

		 double** prop;
		 prop=new double*[n];
		 for(int i=0;i<n;i++){
			 prop[i]=new double[4];
			 for(int j=0;j<4;j++){
				 prop[i][j]=0;
			 }
		 }

		vector<vector<double> > val(4,vector<double>(2));
			 val[0][0]=0;
			 val[0][1]=0;
			 val[1][0]=1;
			 val[1][1]=0;
			 val[2][0]=0;
			 val[2][1]=1;
			 val[3][1]=1;
			 val[3][1]=1;
		for(int i=0;i<n;i++){
			for(int propi=0;propi<4;propi++){
					prop[i][propi]=Test->gis(i,val,propi);
				}
			CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,prop[i][0]+prop[i][1]+prop[i][2]+prop[i][3],0.1);
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

		 GIS *zTest = new GIS(*eX,*esY,*zX,*zY,1,500,0.01,true);
			 //cout << endl;
			 //cout <<zTest->gis(0,0,0,0)-Test->gis(0,0,0,9)  << endl;
			 //cout <<zTest->gis(0,0,1,0)-Test->gis(0,0,1,0) << endl;
			 //cout <<zTest->gis(0,0,0,1)-Test->gis(0,0,0,1) << endl;
			 //cout <<zTest->gis(0,0,1,1)-Test->gis(0,0,1,1) << endl;
			 //cout << endl;
			 //cout <<zTest->gis(1,0,0,0)-Test->gis(1,0,0,5) << endl;
			 //cout <<zTest->gis(1,0,1,0)-Test->gis(1,0,1,0) << endl;
			 //cout <<zTest->gis(1,0,0,1)-Test->gis(1,0,0,0) << endl;
			 //cout <<zTest->gis(1,0,1,1)-Test->gis(1,0,1,1) << endl;
			 //cout << endl;
			 //cout <<zTest->gis(0,1,0,0)-Test->gis(0,1,0,0) << endl;
			 //cout <<zTest->gis(0,1,1,0)-Test->gis(0,1,1,3) << endl;
			 //cout <<zTest->gis(0,1,0,1)-Test->gis(0,1,0,1) << endl;
			 //cout <<zTest->gis(0,1,1,1)-Test->gis(0,1,1,1) << endl;
			 //cout << endl;
			 //cout <<zTest->gis(1,1,0,0)-Test->gis(1,1,0,2) << endl;
			 //cout <<zTest->gis(1,1,1,0)-Test->gis(1,1,1,0) << endl;
			 //cout <<zTest->gis(1,1,0,1)-Test->gis(1,1,0,1) << endl;
			 //cout <<zTest->gis(1,1,1,1)-Test->gis(1,1,1,1) << endl;

			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,0,0)+zTest->gis(0,0,0,1),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,1,0)+zTest->gis(0,0,1,1),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(1,0,0,0)+zTest->gis(1,0,0,1),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(1,0,1,0)+zTest->gis(1,0,1,1),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,1,0,0)+zTest->gis(0,1,0,1),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,1,1,0)+zTest->gis(0,1,1,1),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(1,1,0,0)+zTest->gis(1,1,0,1),0.1);
			 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(1,1,1,0)+zTest->gis(1,1,1,1),0.1);

			 //for(int i=0; i< zTest->getsizeconv()-1; i++ ){
				// cout << zTest->getconv(i)  << " " << zTest->getconv(i+1) << endl;
				// CPPUNIT_ASSERT 	(zTest->getconv(i) >= zTest->getconv(i+1) );
			 //}
			 Test->~GIS();
			 zTest->~GIS();

}
void itTest::SCTwoXTwoY()
{	 srand(time(NULL));
int n=1000;
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

GIS *Test = new GIS(2,*eX,*zX,*zY);

	 Test->setFeatureArraylambda(0,0,1,0,4);
	 Test->setFeatureArraylambda(0,0,1,1,0);
	 Test->setFeatureArraylambda(0,0,0,0,1);
	 Test->setFeatureArraylambda(0,0,0,1,5);

	 Test->setFeatureArraylambda(1,0,1,0,2);
	 Test->setFeatureArraylambda(1,0,1,1,3);
	 Test->setFeatureArraylambda(1,0,0,0,0);
	 Test->setFeatureArraylambda(1,0,0,1,1);

	 Test->setFeatureArraylambda(0,1,1,0,0);
	 Test->setFeatureArraylambda(0,1,1,1,0);
	 Test->setFeatureArraylambda(0,1,0,0,3);
	 Test->setFeatureArraylambda(0,1,0,1,5);

	 Test->setFeatureArraylambda(1,1,1,0,1);
	 Test->setFeatureArraylambda(1,1,1,1,3);
	 Test->setFeatureArraylambda(1,1,0,0,2);
	 Test->setFeatureArraylambda(1,1,0,1,1);

	 double** prop;
	 prop=new double*[n];
	 for(int i=0;i<n;i++){
		 prop[i]=new double[4];
		 for(int j=0;j<4;j++){
			 prop[i][j]=0;
		 }
	 }

	vector<vector<double> > val(4,vector<double>(2));
		 val[0][0]=0;
		 val[0][1]=0;
		 val[1][0]=1;
		 val[1][1]=0;
		 val[2][0]=0;
		 val[2][1]=1;
		 val[3][1]=1;
		 val[3][1]=1;
	for(int i=0;i<n;i++){
		for(int propi=0;propi<4;propi++){
				prop[i][propi]=Test->gis(i,val,propi);
			}
		CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,prop[i][0]+prop[i][1]+prop[i][2]+prop[i][3],0.1);
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

	 SCGIS *zTest = new SCGIS(*eX,*esY,*zX,*zY,1,5000,0.01,true);
		 cout << endl;
		 cout <<zTest->scgis(0,0,0,0)-Test->gis(0,0,0,0)  << endl;
		 cout <<zTest->scgis(0,0,1,0)-Test->gis(0,0,1,0) << endl;
		 cout <<zTest->scgis(0,0,0,1)-Test->gis(0,0,0,1) << endl;
		 cout <<zTest->scgis(0,0,1,1)-Test->gis(0,0,1,1) << endl;
		 cout << endl;
		 cout <<zTest->scgis(1,0,0,0)-Test->gis(1,0,0,0) << endl;
		 cout <<zTest->scgis(1,0,1,0)-Test->gis(1,0,1,0) << endl;
		 cout <<zTest->scgis(1,0,0,1)-Test->gis(1,0,0,1) << endl;
		 cout <<zTest->scgis(1,0,1,1)-Test->gis(1,0,1,1) << endl;
		 cout << endl;
		 cout <<zTest->scgis(0,1,0,0)-Test->gis(0,1,0,0) << endl;
		 cout <<zTest->scgis(0,1,1,0)-Test->gis(0,1,1,0) << endl;
		 cout <<zTest->scgis(0,1,0,1)-Test->gis(0,1,0,1) << endl;
		 cout <<zTest->scgis(0,1,1,1)-Test->gis(0,1,1,1) << endl;
		 cout << endl;
		 cout <<zTest->scgis(1,1,0,0)-Test->gis(1,1,0,0) << endl;
		 cout <<zTest->scgis(1,1,1,0)-Test->gis(1,1,1,0) << endl;
		 cout <<zTest->scgis(1,1,0,1)-Test->gis(1,1,0,1) << endl;
		 cout <<zTest->scgis(1,1,1,1)-Test->gis(1,1,1,1) << endl;

		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(0,0,0,0)+zTest->scgis(0,0,0,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(0,0,1,0)+zTest->scgis(0,0,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(1,0,0,0)+zTest->scgis(1,0,0,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(1,0,1,0)+zTest->scgis(1,0,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(0,1,0,0)+zTest->scgis(0,1,0,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(0,1,1,0)+zTest->scgis(0,1,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(1,1,0,0)+zTest->scgis(1,1,0,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->scgis(1,1,1,0)+zTest->scgis(1,1,1,1),0.1);

		 for(int i=0; i< zTest->getsizeconv()-1; i++ ){
			// cout << zTest->getconv(i)  << " " << zTest->getconv(i+1) << endl;
			 CPPUNIT_ASSERT 	(zTest->getconv(i) >= zTest->getconv(i+1) );
		 }
		 Test->~GIS();
		 zTest->~SCGIS();


}
void itTest::NotBinary(){
	 srand(time(NULL));
		 int n=100;
		 DContainer *eX = new DContainer(n,1);
		 for(int i=0;i< n;i++ ){
			 for(int j=0;j<1;j++){
				 *eX << rand() % 4;
			 }
		 }
		 DContainer *zX = new DContainer(4,1);
		 *zX << 0 << 1 << 2 << 3;
		  DContainer *zY = new DContainer(2,1);
		 *zY << 0 << 1;

		 GIS *Test = new GIS(2,*eX,*zX,*zY);
		 Test->setFeatureArraylambda(0,0,0,0,1);
		 Test->setFeatureArraylambda(0,0,0,1,0);
		 Test->setFeatureArraylambda(0,0,1,0,3);
		 Test->setFeatureArraylambda(0,0,1,1,2);
		 Test->setFeatureArraylambda(0,0,2,0,1);
		 Test->setFeatureArraylambda(0,0,2,1,0);
		 Test->setFeatureArraylambda(0,0,3,0,1);
		 Test->setFeatureArraylambda(0,0,3,1,5);

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
		 GIS *zTest = new GIS(*eX,*esY,*zX,*zY,1,500,0.01,true);

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

		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,0,0)+zTest->gis(0,0,0,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,1,0)+zTest->gis(0,0,1,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,2,0)+zTest->gis(0,0,2,1),0.1);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  1,zTest->gis(0,0,3,0)+zTest->gis(0,0,3,1),0.1);

		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,0,0) ,zTest->gis(0,0,0,0),0.2);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,0,1) ,zTest->gis(0,0,0,1),0.2);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,1,0) ,zTest->gis(0,0,1,0),0.2);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,1,1) ,zTest->gis(0,0,1,1),0.2);

		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,2,0) ,zTest->gis(0,0,2,0),0.2);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,2,1) ,zTest->gis(0,0,2,1),0.2);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,3,0) ,zTest->gis(0,0,3,0),0.2);
		 CPPUNIT_ASSERT_DOUBLES_EQUAL 	(  Test->gis(0,0,3,1) ,zTest->gis(0,0,3,1),0.2);

		 for(int i=0; i< zTest->getsizeconv()-1; i++ ){
			 CPPUNIT_ASSERT 	(zTest->getconv(i)>= zTest->getconv(i+1) );
		 }
		 zTest->~GIS();
		 Test->~GIS();


	}
