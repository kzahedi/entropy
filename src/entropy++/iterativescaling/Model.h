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

        ULContainer* XAlphabet();
        ULContainer* YAlphabet();

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

        double p_y_c_x_d(int yAlphabetIndex, int xAlphabetIndex);
        double p_x_d(int xAlphabetIndex);

        Matrix* p_y_c_x() { return _conditionals;};
        Matrix* p_x()     { return _marginals;};

        int getNrOfUniqueX();
        int getNrOfUniqueY();
        vector<int> getAllColumnsForX();
        vector<int> getAllColumnsForY();

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

        ULContainer* Xdata;
        ULContainer* Ydata;
        ULContainer* Xalphabet;
        ULContainer* Yalphabet;
        ULContainer* _x_alphabet;
        ULContainer* _y_alphabet;

      private:
        int __convertAlphabetToMatrixX(int x);
        int __convertAlphabetToMatrixY(int y);
        ULContainer* __uniqueRows(vector<int> rows, ULContainer* alphabet);

        int           _nrX;
        int           _nrY;

        vector<vector<int> > _Xindices;
        vector<vector<int> > _Yindices;

        // for p(y|x) and p(x)
        vector<int>  _x_indices;
        vector<int>  _y_indices;

        Matrix* _conditionals;
        Matrix* _marginals;
        int     _yAlphabetSize;
    };
  }
}

#endif // __ITERATIVE_SCALING_MODEL_H__
