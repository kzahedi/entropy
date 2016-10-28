#include "feature_test.h"

#include <entropy++/Container.h>
#include <entropy++/MI.h>
#include <entropy++/sparse/MI.h>
#include "../experiments/it/Feature.h"
#include "../experiments/it/GIS.h"
#include "../experiments/it/SCGIS.h"
#include "../experiments/it/Test.h"

#include <iostream>
#include <string>

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( featureTest );
void featureTest::createFeature()
{
	DContainer *aX = new DContainer(3,1);
	(*aX) << 0 << 1 << 2;
	DContainer *aY = new DContainer(2,1);
	(*aY) << 0 << 1;
	//DContainer &aX, DContainer &aY,int colValX, int colValY, int systXsize,int systYsize , double valuelambda
	Feature *test = new Feature(*aX,*aY, 2, 2, 2, 2, -2);
	CPPUNIT_ASSERT(test->getLambdaSize()==0);
	test->setLambda(0,0,3);
	CPPUNIT_ASSERT(test->getLambda(0,0)== 3);
	CPPUNIT_ASSERT(test->getLambda(0,1)== -2);
}
void featureTest::OneFeature()
{
	DContainer *aX = new DContainer(2,1);
	(*aX) << 0 << 1;
	DContainer *aY = new DContainer (2,1);
	(*aY) << 0 << 1;

    vector<vector<int> > systAX(1,vector<int>(0));
    systAX[0].push_back(0);
    vector<vector<int> > systAY(1,vector<int>(0));
    systAY[0].push_back(0);

    IsParameter param;
    param.lambdavalue    = 3.0;
    param.lambdadeltaval = 1.0;
    param.sigma          = 0.01;
    param.maxit          = 500;
    param.konv           = 0.0001;
    param.time           = true;
    param.test           = true;
    param.seconds        = 20;

    DContainer *eX = new DContainer(100,1);
    DContainer *eY = new DContainer(100,1);

    //DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
    GIS *gis = new GIS(*eX, *eY, *aX, *aY, systAX, systAY, param);
    //getFeatureArraylambda(int Feati, int ilambdaX, int ilambdaY)
	CPPUNIT_ASSERT(gis->getFeatureArraylambda(0,0,0) != param.lambdavalue);
	CPPUNIT_ASSERT(gis->getFeatureArraylambda(0,0,1) == param.lambdavalue);
	CPPUNIT_ASSERT(gis->getFeatureArraylambda(0,1,0) == param.lambdavalue);
	CPPUNIT_ASSERT(gis->getFeatureArraylambda(0,1,1) == param.lambdavalue);

	SCGIS *scgis = new SCGIS(*eX, *eY, *aX, *aY, systAX, systAY, param);
	CPPUNIT_ASSERT(scgis->getFeatureArraylambda(0,0,0) != param.lambdavalue);
	CPPUNIT_ASSERT(scgis->getFeatureArraylambda(0,0,1) == param.lambdavalue);
	CPPUNIT_ASSERT(scgis->getFeatureArraylambda(0,1,0) == param.lambdavalue);
	CPPUNIT_ASSERT(scgis->getFeatureArraylambda(0,1,1) == param.lambdavalue);
}
void featureTest::TwoFeatures(){
	DContainer *aX = new DContainer(2,1);
	(*aX) << 0 << 1;
	DContainer *aY = new DContainer (2,1);
	(*aY) << 0 << 1;

    vector<vector<int> > systAX(2,vector<int>(0));
    systAX[0].push_back(0);
    systAX[1].push_back(1);
    vector<vector<int> > systAY(2,vector<int>(0));
    systAY[0].push_back(0);
    systAY[1].push_back(0);

    IsParameter param;
    param.lambdavalue    = 3.0;
    param.lambdadeltaval = 1.0;
    param.sigma          = 0.01;
    param.maxit          = 500;
    param.konv           = 0.0001;
    param.time           = true;
    param.test           = true;
    param.seconds        = 20;
    DContainer *eX = new DContainer(100,2);
    for(int i=0; i<100; i++){
    	(*eX)(i,1)=1;
    }
    DContainer *eY = new DContainer(100,1);
    //DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param)
    GIS *gis = new GIS(*eX, *eY, *aX, *aY, systAX, systAY, param);
	SCGIS *scgis = new SCGIS(*eX, *eY, *aX, *aY, systAX, systAY, param);

		CPPUNIT_ASSERT(gis->getFeatureArraylambda(0,0,0) != param.lambdavalue);
		CPPUNIT_ASSERT(gis->getFeatureArraylambda(0,0,1) == param.lambdavalue);
		CPPUNIT_ASSERT(gis->getFeatureArraylambda(0,1,0) == param.lambdavalue);
		CPPUNIT_ASSERT(gis->getFeatureArraylambda(0,1,1) == param.lambdavalue);
		CPPUNIT_ASSERT(gis->getFeatureArraylambda(1,0,0) == param.lambdavalue);
		CPPUNIT_ASSERT(gis->getFeatureArraylambda(1,0,1) == param.lambdavalue);
		CPPUNIT_ASSERT(gis->getFeatureArraylambda(1,1,0) != param.lambdavalue);
		CPPUNIT_ASSERT(gis->getFeatureArraylambda(1,1,1) == param.lambdavalue);

		CPPUNIT_ASSERT(scgis->getFeatureArraylambda(0,0,0) != param.lambdavalue);
		CPPUNIT_ASSERT(scgis->getFeatureArraylambda(0,0,1) == param.lambdavalue);
		CPPUNIT_ASSERT(scgis->getFeatureArraylambda(0,1,0) == param.lambdavalue);
		CPPUNIT_ASSERT(scgis->getFeatureArraylambda(0,1,1) == param.lambdavalue);
		CPPUNIT_ASSERT(scgis->getFeatureArraylambda(1,0,0) == param.lambdavalue);
		CPPUNIT_ASSERT(scgis->getFeatureArraylambda(1,0,1) == param.lambdavalue);
		CPPUNIT_ASSERT(scgis->getFeatureArraylambda(1,1,0) != param.lambdavalue);
		CPPUNIT_ASSERT(scgis->getFeatureArraylambda(1,1,1) == param.lambdavalue);
}
