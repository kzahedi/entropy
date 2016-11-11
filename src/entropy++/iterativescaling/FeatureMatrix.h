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
        FeatureMatrix(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY,double lambdavalue);
        FeatureMatrix(ULContainer &eX, ULContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY,double lambdavalue);
        FeatureMatrix();
        ~FeatureMatrix();

        vector<int> getMatrixIndexFeat(int i,int j);
        vector<int> getMatrixIndexdX(int i,int j);
        vector<int> getMatrixIndexdY(int i,int j);

      private:
        void __getMatrix(double valuelambda);
        vector<vector<int> > **_mat;
        int  _sizeAlphY;
    };
  }
}

#endif

