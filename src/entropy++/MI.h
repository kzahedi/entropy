#ifndef __MI_H__
#define __MI_H__

#include <entropy++/Container.h>

# define MI_EMPERICAL 2001

class MI
{
  public:

    MI();
    ~MI();

    // MI(X;Y)
    double calculate(Container* X, Container* Y);

  private:

    double __empericalMI(Container* X, Container* Y);

    int _mode;
};


#endif // __MI_H__
