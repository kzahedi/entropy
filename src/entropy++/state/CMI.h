#ifndef __CMI_STATE_H__
#define __CMI_STATE_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace state
  {
    DContainer* CMI(ULContainer* X, ULContainer* Y, ULContainer* Z, int mode);
  }
}


#endif
