#ifndef __CMI_SPARSE_MATRIX_H__
#define __CMI_SPARSE_MATRIX_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace sparse
  {
    double CMI(Container* X, Container* Y, Container* Z, int mode = EMPERICAL);
  }
}

#endif // __CMI_SPARSE_MATRIX_H__
