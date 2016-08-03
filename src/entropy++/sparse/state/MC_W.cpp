#include <entropy++/sparse/state/MC_W.h>

#include <entropy++/sparse/state/CMI.h>

using namespace entropy::sparse::state;

DContainer* entropy::sparse::state::MC_W(DContainer* W2, DContainer* W1, DContainer* A1, int mode)
{
  return CMI(W2, W1, A1, mode);
}
