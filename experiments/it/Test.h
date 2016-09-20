#ifndef __TEST_H__
#define __TEST_H__

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

using namespace std;

class Test
{
  public:
    Test(int colX,int colValY,int rowX, vector<double> lambda,DContainer &aX, DContainer &aY,int maxit, double konv,bool time, int seconds);
    Test(int colX,int colValY,int rowX, vector<double> lambda,DContainer &aX, DContainer &aY);
    // indices and lambda have the same length, indices = feature, lambda = corresponding lambda value
    Test(int colX,int colValY,int rowX, IContainer &indizes, DContainer &lambda, DContainer &aX, DContainer &aY,int maxit, double konv,bool time,int seconds);
    Test(int colX,int colValY,int rowX, DContainer &aX, DContainer &aY,int maxit, double konv,bool time,int seconds,bool comp);
    Test(int colX,int colValY,int rowX, vector<double> lambda,DContainer &aX, DContainer &aY,int maxit, double konv, bool time,bool test,int seconds,int i);
    ~Test();
    vector<double> KL();
    double         KL1();
    void           comparison();
    double         prop(int Feati,int Featj,double ValX,double ValY);
    double         getconv(int ind);
    int            getsizeconv();
    DContainer&    getvalX();
    DContainer&    getvalY();

  private:
    IT                    *_exact;
    GIS                   *_gisTest;
    GISgp                 *_gisgpTest;
    SCGIS                 *_scgisTest;
    SCGISgp               *_scgisgpTest;
    DContainer            *_Y;
    DContainer            *_X;
    DContainer            *_valY;
    DContainer            *_valX;
    int                   _case;
    int                   _sizeColValY;
    int                   _sizeColValX;
    int                   _index;
    vector<double>        _timediff;
    vector<vector<double> > _alphY;
    vector<vector<double> > _alphX;
    bool                  _timetest;

    vector<vector<double> > __getalph(bool valX);

    void      __getValY(int colY,int rowX);
    void      __getValX(int colX,int rowX);
    void      __setlambda(IContainer &indizes, DContainer &values);
    void      __setLambdaRand(vector<double> lambdaval);
    void      __comptime(int maxit, double konv,int seconds);
    void        __fill(vector<double> fill, int i, vector<vector<double> > &Z,bool valX,int rowsAlph,int colval);
};

#endif

