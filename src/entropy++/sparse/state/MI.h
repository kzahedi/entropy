#ifndef __MI_SPARSE_STATE_DEPENDENT_H__
#define __MI_SPARSE_STATE_DEPENDENT_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace sparse
  {
    namespace state
    {
      DContainer* MI(ULContainer* X, ULContainer* Y, int mode = EMPERICAL);
    }
  }
}

#endif
