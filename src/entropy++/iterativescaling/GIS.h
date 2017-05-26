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
          for(vector<Delta*>::const_iterator d = m.deltas.begin(); d != m.deltas.end(); d++)
          {
            cout << **d << endl;
          }
          str << *m._deltaMatcher << endl;
          return str;
        };

      private:
        double         _error;
        vector<double> _s;
        DeltaMatcher*  _deltaMatcher;
    };
  }
}


#endif // __GIS_H__
