#ifndef __MC_MI_STATE_DEPENDENT_H__
#define __MC_MI_STATE_DEPENDENT_H__

#include <entropy++/Container.h>

namespace entropy
{
  namespace state
  {
    Container* MC_MI(Container* W2, Container* W1, Container* S1, Container* A1, int mode = EMPERICAL);
  }
}

#endif
