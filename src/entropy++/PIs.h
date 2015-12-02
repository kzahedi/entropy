#ifndef __PIs_H__
#define __PIs_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

class PIs
{
  public:

    PIs();
    ~PIs();

    // PIs(X;X')
    double calculate(Container* X);

  private:

    double __empericalPIs(Container* X);

    int _mode;
};


#endif // __PIs_H__
