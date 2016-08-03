#include <entropy++/MC_MI.h>

#include <entropy++/MI.h>

double MC_MI(ULContainer* W2, ULContainer* W1, ULContainer* S1, ULContainer* A1, int mode)
{
  double mi_w2_w1 = MI(W2, W1, mode);
  double mi_a1_s1 = MI(S1, A1, mode);
  return mi_w2_w1 - mi_a1_s1;
}
