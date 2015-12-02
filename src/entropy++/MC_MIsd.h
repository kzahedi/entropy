#ifndef __MC_MI_STATE_DEPENDENT_H__
#define __MC_MI_STATE_DEPENDENT_H__

#include <entropy++/Container.h>

class MC_MIsd
{
  public:
    MC_MIsd();

    Container* calculate(Container* W2, Container* W1, Container* S1, Container* A1);

  private:
};

#endif // __MC_MI_H__
