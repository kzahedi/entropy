#include <entropy++/state/MC_W.h>

#include <entropy++/state/CMI.h>

using namespace entropy::state;

Container* entropy::state::MC_W(Container* W2, Container* W1, Container* A1, int mode)
{
  return CMI(W2, W1, A1, mode);
}
