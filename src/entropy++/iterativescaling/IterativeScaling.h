#ifndef __ITERATIVE_SCALING_H__
#define __ITERATIVE_SCALING_H__

#include <entropy++/Csv.h>
#include <entropy++/sparse/MC_MI.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/sparse/state/MC_W.h>
#include <entropy++/sparse/state/MC_MI.h>
#include <entropy++/Matrix.h>
#include <entropy++/SparseMatrix.h>

#include <entropy++/iterativescaling/Feature.h>
#include <entropy++/iterativescaling/ITMatrix.h>
#include <entropy++/iterativescaling/FeatureMatrix.h>
#include <entropy++/iterativescaling/InstanceMatrix.h>
#include <entropy++/iterativescaling/IsParameter.h>

#include <time.h>
#include <string>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

namespace entropy
{
  namespace iterativescaling
  {
    class IterativeScaling
    {

      public:
        IterativeScaling(DContainer &xData,
                         DContainer &yData,
                         DContainer &xAlphabet,
                         DContainer &yAlphabet,
                         vector<vector<int> > systX, // TODO: what does syst mean?
                         vector<vector<int> > systY, 
                         IsParameter param,
                         bool isGis);
        IterativeScaling(ULContainer &xData,
                         ULContainer &yData,
                         DContainer &xAlphabet,
                         DContainer &yAlphabet,
                         vector<vector<int> > systX, // TODO: what does syst mean?
                         vector<vector<int> > systY,
                         IsParameter param,
                         bool isGis);
        IterativeScaling(int ColDataY,
                         DContainer &xData,
                         DContainer &xAlphabet,
                         DContainer &yAlphabet,
                         vector<vector<int> > systX, // TODO: what does syst mean?
                         vector<vector<int> > systY);
        ~IterativeScaling();

        double  prop(int rowX, int rowY);
        double  propAlphX(int indexX, int rowY);
        double  propm(int rowX);
        double  getFeatureArraylambda(int feat, int indexLambdaX, int indexLambdaY);
        void    setFeatureArraylambda(int feat, int indexLambdaX, int indexLambdaY, double valuelambda);
        vector<double> index(int index,bool x, int sizeCol);

      protected:
        double***      __getobs();

        int             _sizeX;
        int             _sizeY;
        int             _sizeColDataX;
        int             _sizeColDataY;
        int             _sizeRowDataX;
        int             _sizeRowDataY;
        int             _sizeSystX;
        double***       _observed;
        bool            _isGis;

        DContainer*     _yAlphabet;
        DContainer*     _xAlphabet;
        DContainer*     _yData;
        DContainer*     _xData;
        FeatureMatrix*  _FM;
        InstanceMatrix* _IM;
        IsParameter     _param;

        vector<vector<int> > _systX;
        vector<vector<int> > _systY;

    };
  }
}

#endif
