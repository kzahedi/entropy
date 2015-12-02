#ifndef __PI_H__
#define __PI_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

class PI
{
  public:

    PI();
    ~PI();

    // PI(X;X')
    double calculate(Container* X);

  private:

    double __empericalPI(Container* X);

    int _mode;
};


#endif // __PI_H__
