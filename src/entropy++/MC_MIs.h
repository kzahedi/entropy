#ifndef __MC_MIs_H__
#define __MC_MIs_H__

#include <entropy++/Container.h>

class MC_MIs
{
  public:
    MC_MIs();

    double calculate(Container* W2, Container* W1, Container* S1, Container* A1);

  private:
};

#endif // __MC_MIs_H__
