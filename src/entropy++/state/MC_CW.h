#ifndef __MC_CW_STATE_DEPENDENT_H__
#define __MC_CW_STATE_DEPENDENT_H__

#include <entropy++/Container.h>

namespace entropy
{
  namespace state
  {
    DContainer* MC_CW(ULContainer* W2, ULContainer* W1, ULContainer* A1, int mode = EMPERICAL);
  }
}

#endif
