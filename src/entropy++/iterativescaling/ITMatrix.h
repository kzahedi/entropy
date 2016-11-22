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
        ITMatrix(ULContainer *eX,
                 ULContainer *eY,
                 ULContainer *aX,
                 ULContainer *aY,
                 ivvector systX,
                 ivvector systY,
                 double lambdavalue);
        virtual ~ITMatrix();
        double  getFeatureArraylambda(int i,int ilambdaX, int ilambdaY);
        double  getFeatureArrayvalue(int i,int rowX, int rowY);
        double  getFeatureArrayvalueAlphY(int feat,int rowX,int indexY);
        double  getFeatureArrayvalueAlphYAlphX(int feat,int indexX,int indexY);
        int 	  getFeatureArraydelta(int i,int indexX, int indexY, int rowDataX, int rowDataY);
        int     getFeatureArraydeltaAlphY(int i,int indexX, int indexY,int rowDataX, int indexDataY);
        int     getFeatureArraydeltaAlphYAlphX(int i,int indexX, int indexY,int indexDataX, int indexDataY);
        void    setFeatureArraylambda(int i, int ilambdaX, int ilambdaY,double valuelambda);
        // vector<int> index(int index,bool x, int sizeCol);
        void index(int* array, int index, bool x, int sizeCol);

      protected:
        void         __featureArray(double valuelambda);

        int          _sizeColDataX;
        int          _sizeColDataY;
        int          _sizeRowDataX;
        int          _sizeRowDataY;
        int          _sizeX;
        int          _sizeY;
        bool         _cmi;
        ivvector     _systX;
        ivvector     _systY;
        ULContainer* _valX;
        ULContainer* _valY;
        ULContainer* _xAlphabet;
        ULContainer* _yAlphabet;
        Feature*     _featureArray;
#if SPEED_OVER_MERMORY

#endif
    };
  }
}

#endif
