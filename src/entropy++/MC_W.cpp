#include <entropy++/MC_W.h>

#include <entropy++/CMI.h>

double MC_W(Container* W2, Container* W1, Container* A1, int mode)
{
  return CMI(W2, W1, A1, mode);
}
