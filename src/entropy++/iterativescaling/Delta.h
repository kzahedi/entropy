#ifndef __DELTA_H__
#define __DELTA_H__

#include <iostream>

namespace entropy
{
  namespace iterativescaling
  {
    class Delta
    {
      public:
        Delta(int xUniqueIndex, int yUniqueIndex);

        bool match(int xUniqueIndex, int yUniqueIndex);

        int getUniqueXIndex();
        int getUniqueYIndex();

        void incObserved();

        void setObserved(double v);

        double observed();

        double expected();

        void setExpected(double v);

        void setLambda(double v);

        double lambda();

        friend std::ostream& operator<<(std::ostream& str, const Delta& m)
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
  }
}

#endif // __DELTA_H__

