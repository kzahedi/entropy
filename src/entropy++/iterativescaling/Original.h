#ifndef __ORIGINAL_H__
#define __ORIGINAL_H__

#include <entropy++/Matrix.h>
#include <entropy++/Container.h>

#include <vector>

using namespace entropy;
using namespace std;

namespace entropy
{
  namespace iterativescaling
  {
    class Original
    {
      public:
        Original();
        Original(int n, vector<vector<int> > features, vector<double> p);
        ~Original();

        void            iterate(int iterations);
       // void   iterate(int klmax);
        double          calculateKL(int iterations);
        vector<double>  getp();
        double          getMarginalProp(int ind,vector<int> feat, vector<double> p);
        double          getConditionalProp(vector<int> featMarg, vector<int> featCond, int ind, vector<double> p);

      private:
        void                _generateAlphabet(int n);
        double              _getprop(vector<double> p, int feat, int ind);
        vector<vector<int> > _features;
        vector<double>      _targetp;
        vector<double>      _p1;
        vector<double>      _p2;
        int                 _sizeAlphabet;
        Matrix*             _alphabet;
    };
  }
}
#endif // __ORIGINAL_H__
