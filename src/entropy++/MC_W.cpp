#include <entropy++/MC_W.h>

#include <entropy++/CMI.h>

double MC_W(ULContainer* W2, ULContainer* W1, ULContainer* A1, int mode)
{
  return CMI(W2, W1, A1, mode);
}
