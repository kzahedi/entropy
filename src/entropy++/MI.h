#ifndef __MI_H__
#define __MI_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  double MI(ULContainer* X, ULContainer* Y, int mode = EMPERICAL);
}

#endif // __MI_H__
