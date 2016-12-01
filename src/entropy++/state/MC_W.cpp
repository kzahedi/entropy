#include <entropy++/state/MC_W.h>

#include <entropy++/state/CMI.h>

using namespace entropy;
using namespace entropy::state;

DContainer* entropy::state::MC_W(ULContainer* W2, ULContainer* W1, ULContainer* A1, int mode)
{
  return CMI(W2, W1, A1, mode);
}
