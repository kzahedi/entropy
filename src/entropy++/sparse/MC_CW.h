#ifndef __MC_CW_SPARSE_MATRIX_H__
#define __MC_CW_SPARSE_MATRIX_H__

#include <entropy++/Container.h>
#include <entropy++/EntropyException.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace sparse
  {
    double MC_CW(ULContainer* W2,
                 ULContainer* W1,
                 ULContainer* A1,
                 int mc_cw_mode = MC_CW_MODE_ENTROPY,
                 int mode       = EMPERICAL)
      throw (EntropyException);
  }
}

#endif
