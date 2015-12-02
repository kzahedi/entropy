#ifndef __CMI_H__
#define __CMI_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

class CMI
{
  public:

    CMI();
    ~CMI();

    // CMI(X;Y|Z)
    double calculate(Container* X, Container* Y, Container* Z);

  private:

    double __empericalCMI(Container* X, Container* Y, Container *Z);

    int _mode;
};


#endif // __CMI_H__
