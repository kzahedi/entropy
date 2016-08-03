#ifndef __MC_W_STATE_H__
#define __MC_W_STATE_H__

#include <entropy++/Entropy.h>
#include <entropy++/Container.h>

namespace entropy
{
  namespace state
  {
    DContainer* MC_W(DContainer* W2, DContainer* W1, DContainer* A1, int mode = EMPERICAL);
  }
}

#endif
