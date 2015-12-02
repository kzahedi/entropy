#ifndef __MI_STATE_DEPENDENT_H__
#define __MI_STATE_DEPENDENT_H__

#include <entropy++/Container.h>

# define MI_EMPERICAL 2001

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
