#ifndef __MC_W_STATE_H__
#define __MC_W_STATE_H__

#include <entropy++/Container.h>

namespace entropy
{
  namespace state
  {
    DContainer* MC_W(ULContainer* W2, ULContainer* W1, ULContainer* A1, int mode = EMPERICAL);
  }
}

#endif
