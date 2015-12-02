#ifndef __MC_MI_STATE_DEPENDENT_SPARSE_MATRIX_H__
#define __MC_MI_STATE_DEPENDENT_SPARSE_MATRIX_H__

#include <entropy++/Container.h>

class MC_MIssd
{
  public:
    MC_MIssd();

    Container* calculate(Container* W2, Container* W1, Container* S1, Container* A1);

  private:
};

#endif
