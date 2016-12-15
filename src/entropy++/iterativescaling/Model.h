#ifndef __ITERATIVE_SCALING_MODEL_H__
#define __ITERATIVE_SCALING_MODEL_H__

#include <entropy++/defs.h>

#include <entropy++/iterativescaling/Feature.h>
#include <entropy++/iterativescaling/Delta.h>
#include <entropy++/Container.h>

#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include <utility> 

using namespace std;

namespace entropy
{
  namespace iterativescaling
  {
    class Model
    {
      public:
        Model();
        ~Model();

        ULContainer* X();
        ULContainer* Y();

        ULContainer* uniqueX(int i);
        ULContainer* uniqueY(int i);

        void createUniqueContainer();

        void setData(ULContainer* X, ULContainer* Y);

        void countObservedFeatures();

        void setFeatures(vector<vector<int> > Xindices,
                         vector<vector<int> > Yindices,
                         vector<Feature*>     features);

        int nrOfFeatures();

        void generateExpected();

        Feature* feature(int i);

        void calculateProbabilities();

        double p_y_c_x(int xUniqueIndex, int yUniqueIndex);
        double p_x(int xUniqueIndex);


      protected:
        vector<Feature*>  features;

      private:

        int           _nrX;
        int           _nrY;

        ULContainer*  _X;
        ULContainer*  _Y;

        ULContainer** _uniqueX;
        ULContainer** _uniqueY;

        vector<vector<int> > _Xindices;
        vector<vector<int> > _Yindices;
    };
  }
}

#endif // __ITERATIVE_SCALING_MODEL_H__
