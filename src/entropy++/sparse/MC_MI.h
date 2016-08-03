#ifndef __MC_MI_SPARSE_MATRIX_H__
#define __MC_MI_SPARSE_MATRIX_H__

#include <entropy++/Container.h>

namespace entropy
{
  namespace sparse
  {
    double MC_MI(ULContainer* W2, ULContainer* W1, ULContainer* S1, ULContainer* A1, int mode = EMPERICAL);
  }
}

#endif
