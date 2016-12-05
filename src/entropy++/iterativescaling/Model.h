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

    typedef std::pair<int,int> relation;

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

        void setRelations(vector<vector<int> > Xindices,
                          vector<vector<int> > Yindices,
                          vector<relation>     relations);

      private:

        int          _nrX;
        int          _nrY;

        ULContainer* _X;
        ULContainer* _Y;

        ULContainer** _uniqueX;
        ULContainer** _uniqueY;

        vector<vector<int> > _Xindices;
        vector<vector<int> > _Yindices;
        vector<relation>     _relations;

    };
  }
}

#endif // __ITERATIVE_SCALING_MODEL_H__
