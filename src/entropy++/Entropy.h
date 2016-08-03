#ifndef __ENTROPY_H__
#define __ENTROPY_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

class Entropy
{
  public:

    Entropy();
    ~Entropy();

    double calulate(IContainer* X);

  private:

    double __empericalEntropy(IContainer* X);

    int _mode;
};


#endif // __ENTROPY_H__
