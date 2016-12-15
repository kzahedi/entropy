#ifndef __KL_H__
#define __KL_H__

#include <entropy++/iterativescaling/Model.h>


namespace entropy
{
  namespace iterativescaling
  {
    class KL
    {
      public:
        KL(Model* p, Model *q);

        double divergence();

      private:
        Model* _p;
        Model* _q;
    };
  }
}

#endif // __KL_H__
