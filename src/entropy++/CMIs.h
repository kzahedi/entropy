#ifndef __CMIs_H__
#define __CMIs_H__

#include <entropy++/Container.h>

# define CMIs_EMPERICAL 2001

class CMIs
{
  public:

    CMIs();
    ~CMIs();

    // CMIs(X;Y|Z)
    double calculate(Container* X, Container* Y, Container* Z);

  private:

    double __empericalCMIs(Container* X, Container* Y, Container *Z);

    int _mode;
};


#endif // __CMIs_H__
