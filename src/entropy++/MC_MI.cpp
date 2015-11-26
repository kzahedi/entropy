#include <entropy++/MC_MI.h>

#include <entropy++/MI.h>

MC_MI::MC_MI()
{ }

double MC_MI::calculate(Container* W2, Container* W1, Container* S1, Container* A1)
{
  MI * mi = new MI();
  double mi_w2_w1 = mi->calculate(W2, W1);
  double mi_a1_s1 = mi->calculate(S1, A1);
  delete mi;
  return mi_w2_w1 = mi_a1_s1;
}
