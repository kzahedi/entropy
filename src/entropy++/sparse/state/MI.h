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
      Container* MI(Container* X, Container* Y, int mode = EMPERICAL);
    }
  }
}

#endif
