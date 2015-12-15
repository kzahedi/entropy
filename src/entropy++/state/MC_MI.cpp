#include <entropy++/state/MC_MI.h>

#include <entropy++/state/MI.h>

using namespace entropy::state;

Container* entropy::state::MC_MI(Container* W2, Container* W1, Container* S1, Container* A1, int mode)
{
  Container *mi_w2_w1 = entropy::state::MI(W2, W1, mode);
  Container *mi_a1_s1 = entropy::state::MI(S1, A1, mode);
  Container *r = new Container(mi_w2_w1->rows(), 1);
  for(int i = 0; i < (int)mi_w2_w1->rows(); i++)
  {
    (*r)(i,0) = mi_w2_w1->get(i,0) - mi_a1_s1->get(i,0);
  }
  return r;
}
