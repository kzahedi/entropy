#ifndef __MC_MI_H__
#define __MC_MI_H__

#include <entropy++/Container.h>

class MC_MI
{
  public:
    MC_MI();

    double calculate(Container* W2, Container* W1, Container* S1, Container* A1);

  private:
};

#endif // __MC_MI_H__
