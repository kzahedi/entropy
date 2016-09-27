#ifndef __CMI_H__
#define __CMI_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  double CMI(ULContainer* X, ULContainer* Y, ULContainer* Z, int mode = EMPERICAL);
}

#endif // __CMI_H__
