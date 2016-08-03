#ifndef __H_SPARSE_MATRIX_H__
#define __H_SPARSE_MATRIX_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace sparse
  {
    double H(DContainer* X, int = EMPERICAL);
  }
}

#endif // __H_SPARSE_MATRIX_H__
