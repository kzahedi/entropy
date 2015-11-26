#include <entropy++/MC_W.h>

#include <entropy++/CMI.h>

MC_W::MC_W()
{
}

double MC_W::calculate(Container* W2, Container* W1, Container* A1)
{
  CMI *cmi = new CMI();
  double r = cmi->calculate(W2, A1, W1);
  delete cmi;
  return r;
}
