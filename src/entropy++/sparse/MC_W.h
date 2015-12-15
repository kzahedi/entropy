#ifndef __MC_W_SPARSE_MATRIX_H__
#define __MC_W_SPARSE_MATRIX_H__

#include <entropy++/Entropy.h>
namespace entropy
{
  namespace sparse
  {
    double MC_W(Container* W2, Container* W1, Container* A1, int mode = EMPERICAL);
  }
}

#endif
