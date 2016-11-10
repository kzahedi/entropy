#ifndef __MC_CW_H__
#define __MC_CW_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <entropy++/EntropyException.h>

namespace entropy
{
  double MC_CW(ULContainer* W2,
               ULContainer* W1,
               ULContainer* A1,
               int mc_cw_mode = MC_CW_MODE_MI,
               int mode       = EMPERICAL)
    throw (EntropyException);
}

#endif // __MC_CW_H__
