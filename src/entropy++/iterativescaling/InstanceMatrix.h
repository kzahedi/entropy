#ifndef __INSTANCEMATRIX_H__
#define __INSTANCEMATRIX_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <entropy++/Matrix.h>
#include <entropy++/iterativescaling/Feature.h>
#include <entropy++/iterativescaling/ITMatrix.h>

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
    class InstanceMatrix : public ITMatrix
    {
      public:
        //InstanceMatrix();
        InstanceMatrix(DContainer  &eX, DContainer  &eY, DContainer &aX, DContainer &aY,ivvector systX, ivvector systY, double valuelambda);
        InstanceMatrix(ULContainer &eX, ULContainer &eY, DContainer &aX, DContainer &aY,ivvector systX, ivvector systY, double valuelambda);
        ~InstanceMatrix();
        // TODO copying vectors can be expensive

        ivector getInstanceMatrixX(int feat, int deltai, int deltaj);
        ivector getInstanceMatrixY(int feat, int deltai, int deltaj);

      private:
        void _getMatrix(double valuelambda);
        // TODO check for a better way to store
        ivvector ***_mat;

    };
  }
}
#endif
