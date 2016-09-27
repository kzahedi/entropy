#ifndef __CONDITIONAL_ENTROPY_SPARSE_MATRIX_H__
#define __CONDITIONAL_ENTROPY_SPARSE_MATRIX_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace sparse
  {
    // H(X|Y)
    double ConditionalEntropy(ULContainer* X, ULContainer* Y, int = EMPERICAL);
  }
}

#endif // __CONDITIONAL_ENTROPY_SPARSE_MATRIX_H__
