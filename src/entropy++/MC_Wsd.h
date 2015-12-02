#ifndef __MC_Wsd_H__
#define __MC_Wsd_H__

#include <entropy++/Entropy.h>
#include <entropy++/Container.h>

class MC_Wsd
{
  public:
    MC_Wsd();

    Container* calculate(Container* W2, Container* W1, Container* A1);

  private:
};

#endif // __MC_Wsd_H__
