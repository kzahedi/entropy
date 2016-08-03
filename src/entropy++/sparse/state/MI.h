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
      DContainer* MI(DContainer* X, DContainer* Y, int mode = EMPERICAL);
    }
  }
}

#endif
