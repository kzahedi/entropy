#include <entropy++/state/MC_MI.h>

#include <entropy++/state/MI.h>

using namespace entropy::state;

DContainer* entropy::state::MC_MI(ULContainer* W2, ULContainer* W1, ULContainer* S1, ULContainer* A1, int mode)
{
  DContainer *mi_w2_w1 = entropy::state::MI(W2, W1, mode);
  DContainer *mi_a1_s1 = entropy::state::MI(S1, A1, mode);
  DContainer *r = new DContainer(mi_w2_w1->rows(), 1);
  for(int i = 0; i < (int)mi_w2_w1->rows(); i++)
  {
    (*r)(i,0) = mi_w2_w1->get(i,0) - mi_a1_s1->get(i,0);
  }
  return r;
}
