#ifndef __MIs_H__
#define __MIs_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

class MIs
{
  public:

    MIs();
    ~MIs();

    // MIs(X;Y)
    double calculate(Container* X, Container* Y);

  private:

    double __empericalMIs(Container* X, Container* Y);

    int _mode;
};


#endif // __MIs_H__
