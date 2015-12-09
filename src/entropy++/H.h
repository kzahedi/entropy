#ifndef __H_H__
#define __H_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

class H
{
  public:

    H();
    ~H();

    // H(X;Y)
    double calculate(Container* X);

  private:

    double __emperical(Container* X);

    int _mode;
};


#endif // __H_H__
