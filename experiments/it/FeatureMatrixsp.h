#ifndef __FEATUREMATRIXSP_H__
#define __FEATUREMATRIXSP_H__

#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include <entropy++/SparseMatrix.h>
#include "Feature.h"

using namespace std;

class FeatureMatrixsp
{
  public:
    FeatureMatrixsp(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue);
    FeatureMatrixsp();
    ~FeatureMatrixsp();
    double  getFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY);
    double  getFeatureArrayvalue(int i, int j,double ValX, double ValY);
    void    setFeatureArraylambda(int i, int j,int ilambdaX, int ilambdaY,double valuelambda);
    int     getFeatureArraydelta(int i, int j,int idelta, int jdelta, double ValX, double ValY);
    int     getMatrixIndexX(int i, int j,int k);
    int     getMatrixIndexY(int i, int j,int k );
    int     getMatrixIndexdX(int i,int j,int k );
    int     getMatrixIndexdY(int i,int j,int k );

  private:
    Feature** FeatureArray(double valuelambda);
    void __getMatrix(double valuelambda);

    int           _sizeColValX;
    int           _sizeColValY;
    int           _sizeRowValX;
    int           _sizeRowValY;
    int           _sizeX;
    int           _sizeY;
    DContainer*   _valX;
    DContainer*   _valY;
    DContainer*   _X;
    DContainer*   _Y;
    Feature**     _FA;
    SparseMatrix* _Feati;
    SparseMatrix* _Featj;
    SparseMatrix* _Delti;
    SparseMatrix* _Deltj;
};

#endif
