#ifndef _GIS_H_
#define _GIS_H_
#include <string>
#include <iostream>
#include <math.h>
#include <vector>
#include "Feature.h"
#include "FeatureMatrix.h"
#include "IT.h"

#include "IsParameter.h"


class GIS : public IT
{
  public:
    GIS(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param);
    ~GIS();
    double getconv(int i);
    int    getsizeconv();
    int    getIterations();

  private:
    double __getFeatconst();
    // get expected frequency of the features, given the current lambdas
    void   __getExpected();
    void   __gis(int maxit, double konv, bool test);
    void   __gis(int maxit, double konv, bool test,int seconds);
    double __calculateIteration(double featconst, bool test);

    double***      _expected;
    double*        _exponent;
    double         _normaliser;
    vector<double> _conv;
    int            _iterations;
};


#endif

