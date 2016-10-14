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
#include "Test.h"
#include "SCGIS.h"
#include "GIS.h"
#include "GISgp.h"
#include "SCGISgp.h"
#include "IT.h"

using namespace std;

class TestMI{
  public:
	TestMI(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, IsParameter param, int version);
	TestMI(ULContainer &eX, ULContainer &eY, int version);
    double            getMI();

  private:
    bool            _cmi;
    DContainer*     _Y;
    DContainer*     _X;
    DContainer*     _valY;
    DContainer*     _valX;
    IT*             _p1;
    IT*             _p2;
    ULContainer*    _valXUL;
    ULContainer*    _valYUL;
};
#endif
