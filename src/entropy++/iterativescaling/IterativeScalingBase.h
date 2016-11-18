#ifndef __ITERATIVE_SCALING_BASE_H__
#define __ITERATIVE_SCALING_BASE_H__

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
    class IterativeScalingBase
    {
      public:
        IterativeScalingBase(ULContainer *xData,
                             ULContainer *yData,
                             ULContainer *xAlphabet,
                             ULContainer *yAlphabet,
                             ivvector systX, // TODO: what does syst mean?
                             ivvector systY,
                             IsParameter param,
                             bool useFeatures);
        IterativeScalingBase(int ColDataY,
                             ULContainer *xData,
                             ULContainer *xAlphabet,
                             ULContainer *yAlphabet,
                             ivvector systX, // TODO: what does syst mean?
                             ivvector systY);
        ~IterativeScalingBase();

        double  prop(int rowX, int rowY);
        double  p_x_c_y(int indexX, int rowY); // TODO: check if correct
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

        ULContainer* _yAlphabet;
        ULContainer* _xAlphabet;
        ULContainer* _yData;
        ULContainer* _xData;
        ITMatrix*   _imatrix;
        IsParameter _param;

        ivvector _systX;
        ivvector _systY;

    };
  }
}

#endif // __ITERATIVE_SCALING_BASE_H__
