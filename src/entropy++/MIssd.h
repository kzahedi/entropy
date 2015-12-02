#ifndef __MIssd_H__
#define __MIssd_H__

#include <entropy++/Container.h>

# define MI_EMPERICAL 2001

class MIssd
{
  public:

    MIssd();
    ~MIssd();

    // MIssd(X;Y)
    Container* calculate(Container* X, Container* Y);

  private:

    Container* __empericalMIssd(Container* X, Container* Y);

    int _mode;
};


#endif // __MIssd_H__
