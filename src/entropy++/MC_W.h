#ifndef __MC_W_H__
#define __MC_W_H__

#include <entropy++/Entropy.h>

class MC_W
{
  public:
    MC_W();
    // ~MC_W();

    //MC_W(const MC_W);
    //MC_W operator=(const MC_W);
    
    // I(W';W|A)
    // W2 = W'
    // W1 = W
    // A1 = A
    double calculate(Container* W2, Container* W1, Container* A1);

  private:
};

#endif // __MC_W_H__
