#ifndef __CMIssd_H__
#define __CMIssd_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

class CMIssd
{
  public:

    CMIssd();
    ~CMIssd();

    // CMIssd(X;Y|Z)
    Container* calculate(Container* X, Container* Y, Container* Z);

  private:

    Container* __emperical(Container* X, Container* Y, Container *Z);

    int _mode;
};


#endif // __CMIssd_H__
