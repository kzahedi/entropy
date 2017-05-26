#ifndef __DELTA_MATCHER_H__
#define __DELTA_MATCHER_H__

#include "Delta.h"
#include <vector>

using namespace std;


namespace entropy
{
  namespace iterativescaling
  {
    typedef vector<Delta*> Dvector;
    class DeltaMatcher
    {
      public:
        DeltaMatcher(int xSize);
        void add(int index, Delta* d);

        vector<Delta*>::iterator begin(int index);
        vector<Delta*>::iterator end(int index);

        friend std::ostream& operator<<(std::ostream& str, const DeltaMatcher& m)
        {
          for(int i = 0; i < m._rows; i++)
          {
            str << "Row " << i << ": " << m._deltas[i].size() << endl;
            for(vector<Delta*>::iterator d = m._deltas[i].begin(); d != m._deltas[i].end(); d++)
            {
              str << "  " << **d << endl;
            }
          }
          return str;
        };

      private:
        Dvector* _deltas;
        int      _rows;
    };
  }
}


#endif // __DELTA_MATCHER_H__
