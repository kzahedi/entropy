#ifndef __DELTA_H__
#define __DELTA_H__

#include <iostream>
#include <vector>

using namespace std;

namespace entropy
{
  namespace iterativescaling
  {
    class Delta
    {
      public:
        Delta(vector<unsigned long> xValues, vector<int> xColumns, vector<unsigned long> yValues, vector<int> yColumns);

        // match piece-wise
        bool matchP(vector<unsigned long>& xValues, vector<unsigned long>& yValues);

        // match full row
        bool matchX(vector<unsigned long>& xValues);

        // match full rows
        bool matchXY(vector<unsigned long>& xValues, vector<unsigned long>& yValues);

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
          str << "X: (";
          for(int i = 0; i < (int)m._xValues.size()-1; i++)
          {
            cout << m._xValues[i] << ",";
          }
          cout << m._xValues[m._xValues.size() - 1] << ") ";
          str << ", Y: (";
          for(int i = 0; i < (int)m._yValues.size()-1; i++)
          {
            cout << m._yValues[i] << ",";
          }
          cout << m._yValues[m._yValues.size() - 1] << ") ";
          str << ", Obs: " << m._observed;
          str << ", Exp: " << m._expected << "]";
          return str;
        };

      private:
        vector<unsigned long> _xValues;
        vector<unsigned long> _yValues;

        vector<int> _xColumns;
        vector<int> _yColumns;

        double _observed;
        double _expected;
        double _lambda;
        double _conditionalProbability;
        double _marginalProbability;
    };
  }
}

#endif // __DELTA_H__

