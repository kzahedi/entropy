#ifndef __ITMATRIX_H__
#define __ITMATRIX_H__


#include <entropy++/defs.h>

#include <entropy++/iterativescaling/Feature.h>

#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

namespace entropy
{
  namespace iterativescaling
  {
    class ITMatrix
    {
      public:
        ITMatrix();
        ITMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY, vector<vector<int> > systX, vector<vector<int> > systY, double lambdavalue);
        ITMatrix(ULContainer &eX, ULContainer &eY, DContainer &aX, DContainer &aY, vector<vector<int> > systX, vector<vector<int> > systY, double lambdavalue);
        virtual ~ITMatrix();
        double  getFeatureArraylambda(int i,int ilambdaX, int ilambdaY);
        double  getFeatureArrayvalue(int i,int rowX, int rowY);
        double  getFeatureArrayvalueAlphY(int feat,int rowX,int indexY);
        double  getFeatureArrayvalueAlphYAlphX(int feat,int indexX,int indexY);
        int 	  getFeatureArraydelta(int i,int indexX, int indexY, int rowValX, int rowValY);
        int     getFeatureArraydeltaAlphY(int i,int indexX, int indexY,int rowValX, int indexValY);
        int     getFeatureArraydeltaAlphYAlphX(int i,int indexX, int indexY,int indexValX, int indexValY);
        void    setFeatureArraylambda(int i, int ilambdaX, int ilambdaY,double valuelambda);
        vector<double> index(int index,bool x, int sizeCol);

      protected:
        void    _FeatureArray(double valuelambda);

        int         _sizeColValX;
        int         _sizeColValY;
        int         _sizeRowValX;
        int         _sizeRowValY;
        int         _sizeX;
        int         _sizeY;
        bool        _cmi;
        vector<vector<int> > _systX;
        vector<vector<int> > _systY;
        DContainer* _valX;
        DContainer* _valY;
        ULContainer* _valXUL;
        ULContainer* _valYUL;
        DContainer* _X;
        DContainer* _Y;
        Feature*    _FA;
    };
  }
}

#endif
