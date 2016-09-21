#ifndef __SCGIS_H__
#define __SCGIS_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <math.h>
#include <entropy++/Matrix.h>
#include "Feature.h"
#include "InstanceMatrix.h"
#include "IT.h"

using namespace std;

class SCGIS : public IT{

public:
  // SCGIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test, bool time,int seconds);
  SCGIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, ItParameter param);
  ~SCGIS();
  double  getconv(int i);
  int     getsizeconv();
  int   getIterations();

private:
  void      __scgis(int maxit, double konv,bool test);
  void      __scgis(int maxit, double konv,bool test,int seconds);
  double    __calculateIteration(bool test);

  int       _iterations;
  double****    _exponent;
  double***   _normaliser;
  double**    _delta;
  vector<double>  _conv;
  ItParameter     _param;

};
#endif
