#ifndef __MC_MI_SPARSE_MATRIX_H__
#define __MC_MI_SPARSE_MATRIX_H__

#include <entropy++/Container.h>

namespace entropy
{
  namespace sparse
  {
    double MC_MI(DContainer* W2, DContainer* W1, DContainer* S1, DContainer* A1, int mode = EMPERICAL);
  }
}

#endif
