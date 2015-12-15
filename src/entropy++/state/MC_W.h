#ifndef __MC_W_STATE_H__
#define __MC_W_STATE_H__

#include <entropy++/Entropy.h>
#include <entropy++/Container.h>

namespace entropy
{
  namespace state
  {
    Container* MC_W(Container* W2, Container* W1, Container* A1, int mode = EMPERICAL);
  }
}

#endif
