#ifndef __FEATURE_H__
#define __FEATURE_H__

#include <entropy++/iterativescaling/Delta.h>

#include <vector>
#include <utility> 

using namespace std;

namespace entropy
{
  namespace iterativescaling
  {

    class Feature : public vector<Delta*>
    {
      public:
        Feature(int xListIndex, int yListIndex);

        int xListIndex();
        int yListIndex();

        void setUniqueXCount(int index, int count);

        int getUniqueXCount(int index);

        void setXAlphabetSize(double x);
        double getXAlphabetSize();

        void setYAlphabetSize(double y);
        double getYAlphabetSize();

        friend std::ostream& operator<<(std::ostream& str, const Feature& m)
        {
          for(vector<Delta*>::const_iterator mf = m.begin(); mf != m.end(); mf++)
          str << **mf << endl;
          return str;
        };

      private:
        int         _xListIndex; // uniqueX index
        int         _yListIndex; // uniqueY index
        double      _xAlphabetSize;
        double      _yAlphabetSize;
        vector<int> _uniqueXCount;
    };
  }
}

#endif // __FEATURE_H__
