#ifndef __GIS_H__
#define __GIS_H__

#include <entropy++/iterativescaling/Model.h>

namespace entropy
{
  namespace iterativescaling
  {
    class GIS : public Model
    { 
      public:
        GIS();

        void   init();
        void   iterate();
        double error();

      private:
        double _error;
    };
  }
}


#endif // __GIS_H__
