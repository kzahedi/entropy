#ifndef _GIS_H_
#define _GIS_H_
#include <string>
#include <iostream>
#include <math.h>
#include <vector>
#include "Feature.h"
#include "FeatureMatrix.h"
#include "IT.h"

class GIS : public IT
{
  public:
    GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue,int maxit, double konv, bool test,bool time,int seconds);
    ~GIS();
    double getconv(int i);
    int    getsizeconv();
    int    getIterations();

  private:
    double __getFeatconst();
    void   __getexp();
    void   __gis(int maxit, double konv, bool test);
    void   __gis(int maxit, double konv, bool test,int seconds);

    double****     _expected;
    double*        _exponent;
    double         _normaliser;
    vector<double> _conv;
    int            _iterations;
};


#endif

