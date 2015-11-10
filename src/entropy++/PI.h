#ifndef __PI_H__
#define __PI_H__

#include <entropy++/Container.h>

# define PI_EMPERICAL 2001

class PI
{
  public:

    PI();
    ~PI();

    // PI(X;X')
    double calulate(Container* X);

  private:

    double __empericalPI(Container* X);

    int _mode;
};


#endif // __PI_H__
