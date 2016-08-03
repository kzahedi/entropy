#include <entropy++/sparse/state/MC_MI.h>

#include <entropy++/sparse/state/MI.h>

using namespace entropy::sparse::state;

DContainer* entropy::sparse::state::MC_MI(DContainer* W2, DContainer* W1, DContainer* S1, DContainer* A1, int mode)
{
  DContainer *mi_w2_w1 = MI(W2, W1, mode);
  DContainer *mi_a1_s1 = MI(S1, A1, mode);
  DContainer *r = new DContainer(mi_w2_w1->rows(), 1);
  for(int i = 0; i < (int)mi_w2_w1->rows(); i++)
  {
    (*r)(i,0) = mi_w2_w1->get(i,0) - mi_a1_s1->get(i,0);
  }
  return r;
}
