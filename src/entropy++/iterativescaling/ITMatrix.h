#ifndef __ITMATRIX_H__
#define __ITMATRIX_H__


#include <entropy++/defs.h>

#include <entropy++/iterativescaling/Feature.h>
#include <entropy++/Container.h>

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
        double  getLambda(int i,int ilambdaX, int ilambdaY);
        double  getValue(int i,int rowX, int rowY);
        double  getValueAlphY(int feat,int rowX,int indexY);
        double  getValueAlphYAlphX(int feat,int indexX,int indexY);
        int     getDelta(int i,int indexX, int indexY, int rowDataX, int rowDataY);
        int     getDeltaAlphY(int i,int indexX, int indexY,int rowDataX, int indexDataY);
        int     getDeltaAlphYAlphX(int i,int indexX, int indexY,int indexDataX, int indexDataY);
        int     getUniqueIndex(int i);
        int     getSizeUnique();
        void    setLambda(int i, int ilambdaX, int ilambdaY,double valuelambda);


#ifdef MEMORY_EFFICIENT
        void index(int* array, int index, bool x, int sizeCol);
#endif


      protected:
        void         __featureArray(double valuelambda);

        int          _sizeColDataX;
        int          _sizeColDataY;
        int          _sizeRowDataX;
        int          _sizeRowDataY;
        int          _sizeX;
        int          _sizeY;
        ivvector     _systX;
        ivvector     _systY;
        ULContainer* _DataX;
        ULContainer* _DataY;
        ULContainer* _xAlphabet;
        ULContainer* _yAlphabet;
        Feature*     _featureArray;
        ULContainer* _UniqueXData;
#ifndef MEMORY_EFFICIENT
        int**  _xFeatureArray; // all possible features for X
        int**  _yFeatureArray; // all possible features for X
        void __fillX();
        void __fillY();
        // double*** _featureArrayvalueAlphY;
#endif
    };
  }
}

#endif
