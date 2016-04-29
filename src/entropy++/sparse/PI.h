#ifndef __PIs_H__
#define __PIs_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace sparse
  {
    double PI(Container* X, int mode = EMPERICAL);
    double PIn(Container* X, int mode = EMPERICAL);
  }
}

#endif // __PIs_H__
