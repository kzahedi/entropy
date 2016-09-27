#ifndef __CONDITIONAL_ENTROPY_H__
#define __CONDITIONAL_ENTROPY_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  double ConditionalEntropy(ULContainer* X, ULContainer* Y, int mode = EMPERICAL);
}


#endif // __CONDITIONAL_ENTROPY_H__
