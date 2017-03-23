#ifndef __SCGIS_H__
#define __SCGIS_H__

#include <entropy++/iterativescaling/IS.h>
#include <entropy++/iterativescaling/RowMatcher.h>

namespace entropy
{
  namespace iterativescaling
  {
    class SCGIS : public IS
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
        RowMatcher* _rowMatcher;
    };
  }
}


#endif // __SCGIS_H__
