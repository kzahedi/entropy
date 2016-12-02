#ifndef __FEATUREMATRIX_H__
#define __FEATUREMATRIX_H__

#include <entropy++/defs.h>
#include <entropy++/Container.h>

#include <entropy++/iterativescaling/Feature.h>
#include <entropy++/iterativescaling/ITMatrix.h>

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
    class FeatureMatrix :public ITMatrix
    {
      public:
        FeatureMatrix(ULContainer *xData,
                      ULContainer *yData,
                      ULContainer *xAlphabet,
                      ULContainer *yAlphabet,
                      ivvector systX,
                      ivvector systY,
                      double lambdavalue);
        FeatureMatrix();
        ~FeatureMatrix();

        // ivector getMatrixIndexFeat(int i,int j);
        // ivector getMatrixIndexdX(int i,int j);
        // ivector getMatrixIndexdY(int i,int j);

        void getMatrixIndexFeat(ivector& r, int i,int j);
        int  getMatrixIndexFeatSize(int i, int j);
        int  getMatrixIndexFeatValue(int i, int j, int k);

        void getMatrixIndexdX(ivector& r, int i, int j);
        int  getMatrixIndexdXValue(int i, int j, int k);
        void getMatrixIndexdY(ivector& r, int i, int j);
        int  getMatrixIndexdYValue(int i, int j, int k);

        int  getUniqueIndex(int i);
        int  getSizeUnique();

      private:
        void __getMatrix(double valuelambda);
        ivvector **  _mat;
        int          _sizeMatrixAlphabetY;
        ULContainer* _UniqueXData;
    };
  }
}

#endif
