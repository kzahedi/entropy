#ifndef __GIS_H__
#define __GIS_H__

#include <entropy++/iterativescaling/IS.h>

namespace entropy
{
  namespace iterativescaling
  {
    class GIS : public IS
    { 
      public:
        GIS();

        void   init();
        void   iterate();
        double error();

      private:
        double         _error;
        vector<double> _s;

    };
  }
}


#endif // __GIS_H__
