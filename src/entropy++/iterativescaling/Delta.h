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
          str << "[";
          str << "X: " << m._xUniqueIndex;
          str << ", Y: " << m._yUniqueIndex;
          str << ", Obs: " << m._observed;
          str << ", Exp: " << m._expected << "]";
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

