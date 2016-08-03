#ifndef __MC_W_STATE_SPARSE_MATRIX_H__
#define __MC_W_STATE_SPARSE_MATRIX_H__

#include <entropy++/Entropy.h>

namespace entropy
{
  namespace sparse
  {
    namespace state
    {
      DContainer* MC_W(DContainer* W2, DContainer* W1, DContainer* A1, int mode = EMPERICAL);
    }
  }
}
#endif
