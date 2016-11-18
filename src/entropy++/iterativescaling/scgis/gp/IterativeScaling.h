#ifndef __SCGISGP_H__
#define __SCGISGP_H__

#include <entropy++/Container.h>
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
      namespace gp
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
            double __calculateIteration(bool test, double sigma);
            void   __scgis(int maxit, double konv,bool test,double sigma);
            void   __scgis(int maxit, double konv,bool test,double sigma,int seconds);

            int             _iterations;
            double***       _exponent;
            double**        _normaliser;
            double***       _delta;
            vector<double>  _conv;
            InstanceMatrix* _im;
        };
      }
    }
  }
}

#endif

