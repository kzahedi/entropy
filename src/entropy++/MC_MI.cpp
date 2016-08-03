#include <entropy++/MC_MI.h>

#include <entropy++/MI.h>

double MC_MI(DContainer* W2, DContainer* W1, DContainer* S1, DContainer* A1, int mode)
{
  double mi_w2_w1 = MI(W2, W1, mode);
  double mi_a1_s1 = MI(S1, A1, mode);
  return mi_w2_w1 - mi_a1_s1;
}
