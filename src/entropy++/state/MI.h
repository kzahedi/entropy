#ifndef __MI_STATE_DEPENDENT_H__
#define __MI_STATE_DEPENDENT_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  namespace state
  {
    Container* MI(Container* X, Container* Y, int mode = EMPERICAL);
  }
}

#endif // __MI_STATE_DEPENDENT_H__
