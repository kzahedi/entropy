#ifndef __ITERATIVE_SCALING_MODEL_H__
#define __ITERATIVE_SCALING_MODEL_H__

#include <entropy++/defs.h>

// #include <entropy++/iterativescaling/Feature.h>
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
    struct MFeature
    {
      public:
        MFeature(int xUniqueIndex, int yUniqueIndex)
        {
          _xUniqueIndex = xUniqueIndex;
          _yUniqueIndex = yUniqueIndex;
          _observed = 0.0;
          _expected = 0.0;
          _lambda   = 1.0;
        }

        bool match(int xUniqueIndex, int yUniqueIndex)
        {
          return (_xUniqueIndex == xUniqueIndex &&
                  _yUniqueIndex == yUniqueIndex);
        };

        int getUniqueXIndex() {return _xUniqueIndex;};
        int getUniqueYIndex() {return _yUniqueIndex;};

        void incObserved()
        {
          _observed += 1.0;
        };

        void setObserved(double v)
        {
          _observed = v;
        }

        double observed()
        {
          return _observed;
        };

        double expected()
        {
          return _expected;
        };

        void setExpected(double v)
        {
          _expected = v;
        }

        void setLambda(double v)
        {
          _lambda = v;
        };

        double lambda()
        {
          return _lambda;
        };

        friend std::ostream& operator<<(std::ostream& str, const MFeature& m)
        {
          str << m._xUniqueIndex << " " << m._yUniqueIndex << " " << m._observed << " " << m._expected;
          return str;
        };

      private:
        int    _xUniqueIndex;
        int    _yUniqueIndex;
        double _observed;
        double _expected;
        double _lambda;
    };


    class Feature : public vector<MFeature*>
    {
      public:
        Feature(int xListIndex, int yListIndex) {_xListIndex = xListIndex; _yListIndex = yListIndex; _alphabetSize = 0.0;};

        int xListIndex() { return _xListIndex;};
        int yListIndex() { return _yListIndex;};

        void setUniqueXCount(int index, int count);

        int getUniqueXCount(int index) {return _uniqueXCount[index];};

        void setRemainingAlphabetSize(double s) {_alphabetSize = s;};
        double getRemainingAlphabetSize() { return _alphabetSize;};

        void setYAlphabetSize(double y){ _yAlphabetSize = y;};
        double yAlphabetSize(){ return _yAlphabetSize;};

        friend std::ostream& operator<<(std::ostream& str, const Feature& m)
        {
          for(vector<MFeature*>::const_iterator mf = m.begin(); mf != m.end(); mf++)
          str << *mf << endl;
          return str;
        };

      private:
        int _xListIndex; // uniqueX index
        int _yListIndex; // uniqueY index
        double _alphabetSize;
        double _yAlphabetSize;
        vector<int> _uniqueXCount;
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

        void setFeatures(vector<vector<int> > Xindices,
                         vector<vector<int> > Yindices,
                         vector<Feature*>     features);

        int nrOfFeatures();

        void generateExpected();

        Feature* feature(int i) {return features[i];};


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
