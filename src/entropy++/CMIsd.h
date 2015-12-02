#ifndef __CMIsd_H__
#define __CMIsd_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

class CMIsd
{
  public:

    CMIsd();
    ~CMIsd();

    Container* calculate(Container* X, Container* Y, Container* Z);

  private:

    Container* __emperical(Container* X, Container* Y, Container *Z);

    int _mode;
};


#endif // __CMIsd_H__
