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

        void   iterate(int iterations);
        void   iterate(double kl);
        double calculateKL(int iterations);
        double calculateKL(vector<double> p, vector<double> q);
        double calculateConditionalKL(vector<double> p, vector<double> q,
                                      vector<int> featMarg, vector<int> featCond);
        double getMarginalProp(int ind,vector<int> feat, vector<double> p);
        double getConditionalProp(vector<int> featMarg, vector<int> featCond,
                                  int ind, vector<double> p);

        vector<double> getp();

        friend std::ostream& operator<<(std::ostream& str, const Original& m)
        {
          int index = 0;
          for(vector<double>::const_iterator p = m._p1.begin(); p != m._p1.end(); p++)
          {
            index++;
            str << index << ": " << *p << endl;
          }
          return str;
        };

      private:
        void                 __generateAlphabet(int n);
        double               __getprop(vector<double>& p, int feat, int ind);
        vector<vector<int> > _features;
        vector<double>       _targetp;
        vector<double>       _p1;
        vector<double>       _p2;
        int                  _sizeAlphabet;
        Matrix*              _alphabet;
    };
  }
}
#endif // __ORIGINAL_H__
