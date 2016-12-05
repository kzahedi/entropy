#ifndef __ITERATIVE_SCALING_MODEL_H__
#define __ITERATIVE_SCALING_MODEL_H__

#include <entropy++/defs.h>

#include <entropy++/iterativescaling/Feature.h>
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

    typedef std::pair<int,int> Relation;

    struct MFeature
    {
      public:
        MFeature(int i, int j, int k, int l)
        {
          _i        = i;
          _j        = j;
          _k        = k;
          _l        = l;
          _observed = 0;
          _lambda   = 0.0;
        }

        bool match(int i, int j, int k, int l)
        {
          return (_i == i && _j == j &&
                  _k == k && _l == l);
        };

        void incObserved()
        {
          _observed++;
        };

        int observed()
        {
          return _observed;
        };

        void setLambda(double v)
        {
          _lambda = v;
        };

        double getLambda()
        {
          return _lambda;
        };

        friend std::ostream& operator<<(std::ostream& str, const MFeature& m)
        {
          str << m._i << " " << m._j << " " << m._k << " " << m._l << " " << m._observed;
          return str;
        };

      private:
        int _i;
        int _j;
        int _k;
        int _l;
        int _observed;
        double _lambda;
    };

    class Model
    {
      public:
        Model();

        ~Model();

        //Model(const Model);
        //Model operator=(const Model);


        ULContainer* X() {return _X;};
        ULContainer* Y() {return _Y;};

        ULContainer* uniqueX(int i) {return _uniqueX[i];};
        ULContainer* uniqueY(int i) {return _uniqueY[i];};

        void createUniqueContainer();

        void setData(ULContainer* X, ULContainer* Y);

        void countObservedFeatures();

        void setRelations(vector<vector<int> > Xindices,
                          vector<vector<int> > Yindices,
                          vector<Relation>     relations);

        int nrOfFeatures();
        MFeature* feature(int i) {return _features[i];};


      private:

        int          _nrX;
        int          _nrY;

        ULContainer* _X;
        ULContainer* _Y;

        ULContainer** _uniqueX;
        ULContainer** _uniqueY;

        vector<vector<int> > _Xindices;
        vector<vector<int> > _Yindices;
        vector<Relation>     _relations;

        vector<MFeature*>    _features;
    };
  }
}

#endif // __ITERATIVE_SCALING_MODEL_H__
