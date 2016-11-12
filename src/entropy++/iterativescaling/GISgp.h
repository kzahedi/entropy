#ifndef __GISGP_H__
#define __GISGP_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <entropy++/Matrix.h>
#include <entropy++/iterativescaling/Feature.h>
#include <entropy++/iterativescaling/InstanceMatrix.h>
#include <entropy++/iterativescaling/GIS.h>
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

    class GISgp : public IterativeScalingBase
    {
      public:
        GISgp(DContainer &xData,
              DContainer &yData,
              DContainer &xAlphabet,
              DContainer &yAlphabet,
              ivvector systX,
              ivvector systY,
              IsParameter param);
        GISgp(ULContainer &xData,
              ULContainer &yData,
              DContainer &xAlphabet,
              DContainer &yAlphabet,
              ivvector systX,
              ivvector systY,
              IsParameter param);
        ~GISgp();

        double  getconv(int i);
        int     getsizeconv();
        int     getIterations();

      private:
        double __calculateIteration(double featconst, double sigma, bool test);
        double __getFeatconst();
        void   __getexp();
        void   __gisgp(int maxit, double konv, double sigma, bool test);
        void   __gisgp(int maxit, double konv, double sigma, bool test, int seconds);

        double***      _expected;
        double***      _delta;
        double**       _exponent;
        double*        _normaliser;
        vector<double> _conv;
        int            _iterations;
        FeatureMatrix* _fm;
    };
  }
}

#endif

