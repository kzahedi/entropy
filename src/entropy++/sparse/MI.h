#ifndef __MIs_H__
#define __MIs_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace sparse
  {
    double MI(ULContainer* X, ULContainer* Y, int mode = EMPERICAL);
  }
}

#endif // __MIs_H__
