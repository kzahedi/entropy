#ifndef __CMI_STATE_SPARSE_MATRIX_H__
#define __CMI_STATE_SPARSE_MATRIX_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace sparse
  {
    namespace state
    {
      Container* CMI(Container* X, Container* Y, Container* Z, int mode = EMPERICAL);
    }
  }
}


#endif 
