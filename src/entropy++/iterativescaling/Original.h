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
        ~Original();

        void setData(ULContainer* X);
        void addFeature(vector<int> f);
        void init();
        void iterate();

      private:
        Matrix* _joint;       // (n+1)
        Matrix* _conditional; // (n)
        Matrix* _emperical;   // from data

        ULContainer* _X;
        ULContainer* _Xalphabet;
        Matrix** _marginal; // marginal for each feature


        vector<vector<int> > _features;
        int  _n;
        bool _initCalled;

    };
  }
}
#endif // __ORIGINAL_H__
