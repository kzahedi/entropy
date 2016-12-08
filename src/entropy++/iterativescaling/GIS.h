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
        // ~GIS : public Model();

        //GIS : public Model(const GIS : public Model);
        //GIS : public Model operator=(const GIS : public Model);

        void init();
        void iterate();
        double error();

      private:
        double _error;
    };
  }
}


#endif // __GIS_H__
