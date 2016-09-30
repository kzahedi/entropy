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
    double            getMI();

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
