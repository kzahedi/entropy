
#ifndef __Joint_GIS_H__
#define __Joint_GIS_H__

#include <entropy++/iterativescaling/IS.h>
#include <entropy++/iterativescaling/DeltaMatcher.h>

namespace entropy
{
  namespace iterativescaling
  {
    class Joint_GIS : public IS
    {
      public:
        Joint_GIS();
        ~Joint_GIS();
        void            setFeatures(vector<vector<int> > Xindices, vector<Feature*> f);
        void            countObservedFeatures();
        void            init();
        void            setAlphabet(ULContainer* XAlphabet);
        void            iterate();
        double          error();
        void            calculateProbabilities();
        bool            matchX(vector<unsigned long> x, vector<unsigned long> y);
        vector<double>  getProbabilities();
        double          calculateKL(vector<double> q);

        friend std::ostream& operator<<(std::ostream& str, const Joint_GIS& m)
        {
          cout << "Data:" << endl;
          for(int i = 0; i < m.Xdata->rows(); i++)
          {
            cout << m.Xdata->get(i,0);
            for(int j = 1; j < m.Xdata->columns(); j++)
            {
              cout << "," << m.Xdata->get(i,j);
            }
            cout << endl;
          }
          for(vector<Delta*>::const_iterator d = m.deltas.begin(); d != m.deltas.end(); d++)
          {
            cout << **d << endl;
          }
          str << *m._deltaMatcher << endl;
          return str;
        };

      private:
        double                _error;
        double                _sizeAlphabet;
        vector<double>        _marginals;
        vector<double>        _s;
        vector<int>           _indices;
        DeltaMatcher*         _deltaMatcher;
        vector<vector<int> >  _Indices;
    };
  }
}


#endif
