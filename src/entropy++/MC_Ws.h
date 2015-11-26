#ifndef __MC_Ws_H__
#define __MC_Ws_H__

#include <entropy++/Entropy.h>

class MC_Ws
{
  public:
    MC_Ws();
    // ~MC_Ws();

    //MC_Ws(const MC_Ws);
    //MC_Ws operator=(const MC_Ws);
    
    // I(W';W|A)
    // W2 = W'
    // W1 = W
    // A1 = A
    double calculate(Container* W2, Container* W1, Container* A1);

  private:
};

#endif // __MC_Ws_H__
