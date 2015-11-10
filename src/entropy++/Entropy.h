#ifndef __ENTROPY_H__
#define __ENTROPY_H__

#include <entropy++/Container.h>

# define Entropy_EMPERICAL 2001

class Entropy
{
  public:

    Entropy();
    ~Entropy();

    double calulate(Container* X);

  private:

    double __empericalEntropy(Container* X);

    int _mode;
};


#endif // __ENTROPY_H__
