#include <entropy++/MC_MIssd.h>

#include <entropy++/MIssd.h>

MC_MIssd::MC_MIssd()
{ }

Container* MC_MIssd::calculate(Container* W2, Container* W1, Container* S1, Container* A1)
{
  MIssd* mi = new MIssd();
  Container *mi_w2_w1 = mi->calculate(W2, W1);
  Container *mi_a1_s1 = mi->calculate(S1, A1);
  Container *r = new Container(mi_w2_w1->rows(), 1);
  for(int i = 0; i < (int)mi_w2_w1->rows(); i++)
  {
    (*r)(i,0) = mi_w2_w1->get(i,0) - mi_a1_s1->get(i,0);
  }
  return r;
}
