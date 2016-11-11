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
                         ivvector systX, // TODO: what does syst mean?
                         ivvector systY, 
                         IsParameter param,
                         bool useFeatures);
        IterativeScaling(ULContainer &xData,
                         ULContainer &yData,
                         DContainer &xAlphabet,
                         DContainer &yAlphabet,
                         ivvector systX, // TODO: what does syst mean?
                         ivvector systY,
                         IsParameter param,
                         bool useFeatures);
        IterativeScaling(int ColDataY,
                         DContainer &xData,
                         DContainer &xAlphabet,
                         DContainer &yAlphabet,
                         ivvector systX, // TODO: what does syst mean?
                         ivvector systY);
        ~IterativeScaling();

        double  prop(int rowX, int rowY);
        double  propAlphX(int indexX, int rowY);
        double  propm(int rowX);
        double  getFeatureArraylambda(int feat, int indexLambdaX, int indexLambdaY);
        void    setFeatureArraylambda(int feat, int indexLambdaX, int indexLambdaY, double valuelambda);
        dvector index(int index,bool x, int sizeCol);

      protected:
        double***   __getobs();

        int         _sizeX;
        int         _sizeY;
        int         _sizeColDataX;
        int         _sizeColDataY;
        int         _sizeRowDataX;
        int         _sizeRowDataY;
        int         _sizeSystX;
        double***   _observed;

        DContainer* _yAlphabet;
        DContainer* _xAlphabet;
        DContainer* _yData;
        DContainer* _xData;
        ITMatrix*   _im;
        IsParameter _param;

        ivvector _systX;
        ivvector _systY;

    };
  }
}

#endif
