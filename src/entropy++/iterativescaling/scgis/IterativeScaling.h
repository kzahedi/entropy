#ifndef __IterativeScaling_H__
#define __IterativeScaling_H__

#include <entropy++/defs.h>
#include <entropy++/Matrix.h>
#include <entropy++/iterativescaling/Feature.h>
#include <entropy++/iterativescaling/InstanceMatrix.h>
#include <entropy++/iterativescaling/IterativeScalingBase.h>

#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <math.h>

using namespace std;

namespace entropy
{
  namespace iterativescaling
  {
    namespace scgis
    {
      class IterativeScaling : public IterativeScalingBase
      {
        public:
          IterativeScaling(ULContainer *xData,
                ULContainer *yData,
                ULContainer *xAlphabet,
                ULContainer *yAlphabet,
                ivvector systX,
                ivvector systY,
                IsParameter param);
          ~IterativeScaling();
          double  getconv(int i);
          int     getsizeconv();
          int     getIterations();

        private:
          void           __scgis(int maxit, double konv, bool test);
          void           __scgis(int seconds, bool test);
          void           __scgis(double konv,int seconds, bool test);
          double         __calculateIteration(bool test);

          double         _delta;
          int            _iterations;
          double***      _exponent;
          double**       _normaliser;

          vector<double> _conv;
          IsParameter    _param;
          InstanceMatrix* _im;
      };
    }
  }
}

#endif
