#ifndef __SCGIS_H__
#define __SCGIS_H__

#include <entropy++/iterativescaling/Model.h>

namespace entropy
{
  namespace iterativescaling
  {
    class SCGIS : public Model
    { 
      public:
        SCGIS();
        ~SCGIS();

        void   init();
        void   iterate();
        double error();

      private:
        double _error;
        Matrix* _z;
        Matrix* _s;
    };
  }
}


#endif // __SCGIS_H__
