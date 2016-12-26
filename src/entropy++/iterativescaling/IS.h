#ifndef __IS_H__
#define __IS_H__

#include <entropy++/iterativescaling/Model.h>

namespace entropy
{
  namespace iterativescaling
  {
    class IS : public Model
    { 
      public:
        IS() : Model() {};
        virtual void   init()    = 0;
        virtual void   iterate() = 0;
        virtual double error()   = 0;
    };
  }
}


#endif // __IS_H__
