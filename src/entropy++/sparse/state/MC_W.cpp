#include <entropy++/sparse/state/MC_W.h>

#include <entropy++/sparse/state/CMI.h>

using namespace entropy::sparse::state;

DContainer* entropy::sparse::state::MC_W(ULContainer* W2, ULContainer* W1, ULContainer* A1, int mode)
{
  return CMI(W2, W1, A1, mode);
}
