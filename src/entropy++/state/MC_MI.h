#ifndef __MC_MI_STATE_DEPENDENT_H__
#define __MC_MI_STATE_DEPENDENT_H__

#include <entropy++/Container.h>

namespace entropy
{
  namespace state
  {
    DContainer* MC_MI(DContainer* W2, DContainer* W1, DContainer* S1, DContainer* A1, int mode = EMPERICAL);
  }
}

#endif
