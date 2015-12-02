#include <entropy++/MC_Wssd.h>

#include <entropy++/CMIssd.h>

MC_Wssd::MC_Wssd()
{ }

Container* MC_Wssd::calculate(Container* W2, Container* W1, Container* A1)
{
  CMIssd*    cmi = new CMIssd();
  Container* r   = cmi->calculate(W2, W1, A1);
  delete cmi;
  return r;
}
