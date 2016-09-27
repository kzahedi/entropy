#include <entropy++/MC_CW.h>

#include <entropy++/MI.h>

double MC_CW(ULContainer* W2, ULContainer* W1, ULContainer* A1, int mode)
{
  double mi_w2_w1 = MI(W2, W1, mode);
  double mi_w2_a1 = MI(W2, A1, mode);
  return mi_w2_w1 - mi_w2_a1;
}
