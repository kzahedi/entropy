#ifndef __TESTMI_H__
#define __TESTMI_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <math.h>
#include <entropy++/Matrix.h>
#include "InstanceMatrix.h"
#include "SCGIS.h"
#include "GIS.h"
#include "GISgp.h"
#include "SCGISgp.h"
#include "IT.h"

using namespace std;

class TestMI{
	//int colX,int colValY, int rowX, vector<double> lambda, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param);
  public:
	TestMI(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<double> lambda, vector<vector<int> > systAX, vector<vector<int> > systAY, vector<vector<int> > systBX, vector<vector<int> > systBY, IsParameter param, int version);
    TestMI(DContainer &eX, DContainer &eY, IContainer &indizes, DContainer &lambdas, DContainer &aX, DContainer &aY, vector<vector<int> > systAX, vector<vector<int> > systAY, vector<vector<int> > systBX, vector<vector<int> > systBY, IsParameter param, int version);
    double getMI();

  private:
    DContainer*     _Y;
    DContainer*     _X;
    DContainer*     _valY;
    DContainer*     _valX;
    GIS*            _p1;
    GIS*            _p2;
};
#endif
