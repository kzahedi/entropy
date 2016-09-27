#ifndef __MC_MI_H__
#define __MC_MI_H__

#include <entropy++/Container.h>

namespace entropy
{
  double MC_MI(ULContainer* W2, ULContainer* W1, ULContainer* S1, ULContainer* A1, int mode = EMPERICAL);
}

#endif // __MC_MI_H__
