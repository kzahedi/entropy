#ifndef __GIS_H__
#define __GIS_H__

#include <entropy++/iterativescaling/IS.h>
#include <entropy++/iterativescaling/DeltaMatcher.h>

namespace entropy
{
  namespace iterativescaling
  {
    class GIS : public IS
    {
      public:
        GIS();
        ~GIS();
        void   init();
        void   iterate();
        double error();

        friend std::ostream& operator<<(std::ostream& str, const GIS& m)
        {
          cout << "Data:" << endl;
          for(int i = 0; i < m.Xdata->rows(); i++)
          {
            cout << m.Xdata->get(i,0);
            for(int j = 1; j < m.Xdata->columns(); j++)
            {
              cout << "," << m.Xdata->get(i,j);
            }
            cout << " : ";
            cout << m.Ydata->get(i,0);
            for(int j = 1; j < m.Ydata->columns(); j++)
            {
              cout << "," << m.Ydata->get(i,j);
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
        void           __unsetDeltas();
        double         _error;
        vector<double> _s;
        DeltaMatcher*  _deltaMatcher;
    };
  }
}


#endif // __GIS_H__
