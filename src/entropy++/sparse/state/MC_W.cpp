#include <entropy++/sparse/state/MC_W.h>

#include <entropy++/sparse/state/CMI.h>

using namespace entropy::sparse::state;

Container* entropy::sparse::state::MC_W(Container* W2, Container* W1, Container* A1, int mode)
{
  return CMI(W2, W1, A1, mode);
}
