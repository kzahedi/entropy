#ifndef __CMI_STATE_H__
#define __CMI_STATE_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace state
  {
    DContainer* CMI(DContainer* X, DContainer* Y, DContainer* Z, int mode);
  }
}


#endif
