#include <entropy++/MC_Wsd.h>

#include <entropy++/CMIsd.h>

MC_Wsd::MC_Wsd()
{ }

Container* MC_Wsd::calculate(Container* W2, Container* W1, Container* A1)
{
  CMIsd*     cmi = new CMIsd();
  Container* r   = cmi->calculate(W2, W1, A1);
  delete cmi;
  return r;
}
