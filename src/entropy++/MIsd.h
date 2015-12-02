#ifndef __MI_STATE_DEPENDENT_H__
#define __MI_STATE_DEPENDENT_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

class MIsd
{
  public:

    MIsd();
    ~MIsd();

    // MI(X;Y)
    Container* calculate(Container* X, Container* Y);

  private:

    Container* __empericalMI(Container* X, Container* Y);

    int _mode;
};


#endif // __MI_STATE_DEPENDENT_H__
