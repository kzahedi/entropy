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

        void setExpected(double v);
        double expected();

        void setLambda(double v);
        double lambda();

        void setConditionalProbability(double value);
        double conditionalProbability();

        void setMarginalProbability(double value);
        double marginalProbability();

        friend std::ostream& operator<<(std::ostream& str, const Delta& m)
        {
          str << "Delta:" << std::endl;
          str << "  X-index:  " << m._xUniqueIndex << std::endl;
          str << "  Y-index:  " << m._yUniqueIndex << std::endl;
          str << "  Observed: " << m._observed << std::endl;
          str << "  Expected: " << m._expected << std::endl;
          return str;
        };

      private:
        int    _xUniqueIndex;
        int    _yUniqueIndex;
        double _observed;
        double _expected;
        double _lambda;
        double _conditionalProbability;
        double _marginalProbability;
    };
  }
}

#endif // __DELTA_H__

