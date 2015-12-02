#include <entropy++/MC_Ws.h>

#include <entropy++/CMIs.h>

MC_Ws::MC_Ws()
{
}

double MC_Ws::calculate(Container* W2, Container* W1, Container* A1)
{
  CMIs *cmi = new CMIs();
  double r = cmi->calculate(W2, W1, A1);
  delete cmi;
  return r;
}
