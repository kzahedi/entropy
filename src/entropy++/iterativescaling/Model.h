#ifndef __ITERATIVE_SCALING_MODEL_H__
#define __ITERATIVE_SCALING_MODEL_H__

#include <entropy++/defs.h>

#include <entropy++/iterativescaling/Feature.h>
#include <entropy++/iterativescaling/Delta.h>
#include <entropy++/Container.h>
#include <entropy++/Matrix.h>

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

        void createUniqueContainer();

        void setData(ULContainer* X, ULContainer* Y);

        void countObservedFeatures();

        void setFeatures(vector<vector<int> > Xindices,
                         vector<vector<int> > Yindices,
                         vector<Feature*>     features);

        int nrOfFeatures();

        Feature* feature(int i);

        void calculateProbabilities();

        double p_y_c_x(int yUniqueIndex, int xUniqueIndex);
        double p_x(int xUniqueIndex);

        int getNrOfUniqueX();
        int getNrOfUniqueY();

        friend std::ostream& operator<<(std::ostream& str, const Model& m)
        {
          int index = 0;
          for(vector<Delta*>::const_iterator d = m.deltas.begin(); d != m.deltas.end(); d++)
          {
            index++;
            str << "Delta " << index << ": " << **d << endl;
          }
          return str;
        };

      protected:
        vector<Feature*> features;
        vector<Delta*>   deltas;

        ULContainer*  Xdata;
        ULContainer*  Ydata;
        ULContainer*  Xalphabet;
        ULContainer*  Yalphabet;

      private:

        int           _nrX;
        int           _nrY;

        vector<vector<int> > _Xindices;
        vector<vector<int> > _Yindices;

        Matrix* _conditionals;
        Matrix* _marginals;
        int     _yAlphabetSize;
    };
  }
}

#endif // __ITERATIVE_SCALING_MODEL_H__
