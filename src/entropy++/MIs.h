#ifndef __MIs_H__
#define __MIs_H__

#include <entropy++/Container.h>

# define MIs_EMPERICAL 2001

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
