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
  public:
	TestMI(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, vector<vector<int> > systAX, vector<vector<int> > systAY, vector<vector<int> > systBX, vector<vector<int> > systBY, IsParameter param, int version);
    TestMI(int colX,int colValY, int rowX,DContainer &aX, DContainer &aY,vector<double> lambda, vector<vector<int> > systAX, vector<vector<int> > systAY, vector<vector<int> > systBX, vector<vector<int> > systBY,vector<vector<int> > systCX, vector<vector<int> > systCY, IsParameter param);
    double getMI();

  private:
    DContainer*     _Y;
    DContainer*     _X;
    DContainer*     _valY;
    DContainer*     _valX;
    IT*             _p1;
    IT*             _p2;
    IT*	            _exact;
};
#endif
