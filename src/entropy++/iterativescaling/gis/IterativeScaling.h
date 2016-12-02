#ifndef _IterativeScaling_H_
#define _IterativeScaling_H_

#include <entropy++/iterativescaling/Feature.h>
#include <entropy++/iterativescaling/FeatureMatrix.h>
#include <entropy++/iterativescaling/IterativeScalingBase.h>
#include <entropy++/iterativescaling/IsParameter.h>
#include <entropy++/Container.h>

#include <string>
#include <iostream>
#include <math.h>
#include <vector>


namespace entropy
{
  namespace iterativescaling
  {
    namespace gis
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

          double getconv(int i);
          int    getsizeconv();
          int    getIterations();

        private:
          double __getFeatconst();
          // get expected frequency of the features, given the current lambdas
          void   __getExpected();
          void   __gis(int maxit, double konv, bool test);
          void   __gis(int seconds, bool test);
          void   __gis(double konv, int seconds, bool test);
          double __calculateIteration(double featconst, bool test);

          double***      _expected;
          double**       _exponent;
          double*        _normaliser;
          int            _iterations;
          dvector        _conv;
          FeatureMatrix* _fm;
      };
    }
  }
}


#endif

