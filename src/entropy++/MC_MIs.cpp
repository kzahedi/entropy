#include <entropy++/MC_MIs.h>

#include <entropy++/MIs.h>

MC_MIs::MC_MIs()
{ }

double MC_MIs::calculate(Container* W2, Container* W1, Container* S1, Container* A1)
{
  MIs * mi = new MIs();
  double mi_w2_w1 = mi->calculate(W2, W1);
  double mi_a1_s1 = mi->calculate(S1, A1);
  delete mi;
  return mi_w2_w1 - mi_a1_s1;
}
