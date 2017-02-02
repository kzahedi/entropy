#ifndef __H_H__
#define __H_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  double H(ULContainer* X, int mode = EMPERICAL);
}

#endif // __H_H__
