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
    Test(int colX,int colValY, int rowX, vector<double> lambda, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param);
    Test(int colX,int colValY, int rowX, vector<double> lambda,DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY);
    Test(int colX,int colValY, int rowX, IContainer &indizes, DContainer &lambda ,DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param);
    Test(int colX,int colValY, int rowX, vector<double> lambda, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param, int i);
    ~Test();
    vector<double>  KL();
    double          KL1();
    void            comparison();
    double          prop( int indexX, int indexY);
    double          getconv(int ind);
    int             getsizeconv();
    vector<double>  propAll(int indexX, int indexY);
    vector<double>  getconvAll(int ind);
    vector<int>     getsizeconvAll();
    DContainer&     getvalX();
    DContainer&     getvalY();

  private:
    IT*                     _exact;
    GIS*                    _gisTest;
    GISgp*                  _gisgpTest;
    SCGIS*                  _scgisTest;
    SCGISgp*                _scgisgpTest;
    DContainer*             _Y;
    DContainer*             _X;
    DContainer*             _valY;
    DContainer*             _valX;
    int                     _case;
    int                     _sizeColValY;
    int                     _sizeColValX;
    int                     _index;
    vector<double>          _timediff;
    vector<vector<int> >    _systX;
    vector<vector<int> >    _systY;
    bool                    _timetest;

    void __getValY(int colY,int rowX);
    void __getValX(int colX,int rowX);
    void __setlambda(IContainer &indizes, DContainer &values);
    void __setLambdaRand(vector<double> lambdaval);
    void __comptime(IsParameter param);
};
#endif

