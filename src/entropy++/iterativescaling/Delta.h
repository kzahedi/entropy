#ifndef __DELTA_H__
#define __DELTA_H__

#include <iostream>
#include <vector>
#include <boost/thread.hpp>

using namespace std;

namespace entropy
{
  namespace iterativescaling
  {
    class Delta
    {
      public:
        Delta(vector<unsigned long> xValues, vector<int> xIndices, vector<unsigned long> yValues, vector<int> yIndices);

        // match piece-wise
        bool match(vector<unsigned long> xValues, vector<unsigned long> yValues);

        // match full row
        bool matchX(vector<unsigned long> xValues);

        bool matchXY(vector<unsigned long> xValues, vector<unsigned long> yValues);

        bool matchY(vector<unsigned long> yValues);

        void incObserved();

        void setObserved(double v);
        double observed();

        void updateExpected(double v);
        void setExpected(double v);
        double expected();

        void updateLambda(double v);
        void setLambda(double v);
        double lambda();

        void setConditionalProbability(double value);
        double conditionalProbability();

        void setMarginalProbability(double value);
        double marginalProbability();

        void setInputOnly(); // if used, then this delta will only match X
        bool isInputOnly();

        void setOutputOnly(); // if used, then this delta will only match X
        bool isOutputOnly();

        friend std::ostream& operator<<(std::ostream& str, const Delta& m)
        {
          str << "X: (";
          for(int i = 0; i < (int)m._xValues.size()-1; i++)
          {
            str << m._xValues[i] << ",";
          }
          if(m._xValues.size() > 0) str << m._xValues[m._xValues.size() - 1];
          str << ") [";
          for(int i = 0; i < (int)m._xIndices.size()-1; i++)
          {
            str << m._xIndices[i] << ",";
          }
          if(m._xIndices.size() > 0) str << m._xIndices[m._xIndices.size() - 1];
          str << "], Y: (";
          for(int i = 0; i < (int)m._yValues.size()-1; i++)
          {
            str << m._yValues[i] << ",";
          }
          if(m._yValues.size() > 0) str << m._yValues[m._yValues.size() - 1];
          str << ") [";
          for(int i = 0; i < (int)m._yIndices.size()-1; i++)
          {
            str << m._yIndices[i] << ",";
          }
          if(m._yIndices.size() > 0) str << m._yIndices[m._yIndices.size() - 1];
          str << "]";
          str << ", Obs: " << m._observed;
          str << ", Exp: " << m._expected;
          str << ", Lambda: " << m._lambda;
          return str;
        };

      private:
        vector<unsigned long> _xValues;
        vector<unsigned long> _yValues;

        vector<int> _xIndices;
        vector<int> _yIndices;

        double       _observed;
        double       _expected;
        double       _lambda;
        double       _conditionalProbability;
        double       _marginalProbability;
        bool         _inputOnly;
        bool         _outputOnly;
        boost::mutex _mutex;
    };
  }
}

#endif // __DELTA_H__

