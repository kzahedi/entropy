#ifndef __Hs_H__
#define __Hs_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

class Hs
{
  public:

    Hs();
    ~Hs();

    // Hs(X;Y)
    double calculate(Container* X);

  private:

    double __emperical(Container* X);

    int _mode;
};


#endif // __Hs_H__
