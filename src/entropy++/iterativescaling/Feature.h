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

        friend std::ostream& operator<<(std::ostream& str, const Feature& m)
        {
          for(vector<Delta*>::const_iterator d = m.begin(); d != m.end(); d++) str << **d << endl;
          return str;
        };

      private:
        int         _xListIndex; // uniqueX index
        int         _yListIndex; // uniqueY index
    };
  }
}

#endif // __FEATURE_H__
