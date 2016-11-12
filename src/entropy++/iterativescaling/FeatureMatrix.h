#ifndef __FEATUREMATRIX_H__
#define __FEATUREMATRIX_H__

#include <entropy++/defs.h>

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
        FeatureMatrix(DContainer *xData,
                      DContainer *yData,
                      DContainer *xAlphabet,
                      DContainer *yAlphabet,
                      ivvector systX,
                      ivvector systY,
                      double lambdavalue);
        FeatureMatrix(ULContainer *xData,
                      ULContainer *yData,
                      DContainer *xAlphabet,
                      DContainer *yAlphabet,
                      ivvector systX,
                      ivvector systY,
                      double lambdavalue);
        FeatureMatrix();
        ~FeatureMatrix();

        ivector getMatrixIndexFeat(int i,int j);
        ivector getMatrixIndexdX(int i,int j);
        ivector getMatrixIndexdY(int i,int j);

      private:
        void __getMatrix(double valuelambda);
        ivvector **_mat;
        int  _sizeAlphY;
    };
  }
}

#endif
