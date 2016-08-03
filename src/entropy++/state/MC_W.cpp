#include <entropy++/state/MC_W.h>

#include <entropy++/state/CMI.h>

using namespace entropy::state;

DContainer* entropy::state::MC_W(DContainer* W2, DContainer* W1, DContainer* A1, int mode)
{
  return CMI(W2, W1, A1, mode);
}
