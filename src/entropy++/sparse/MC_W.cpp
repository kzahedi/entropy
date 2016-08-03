#include <entropy++/sparse/MC_W.h>

#include <entropy++/sparse/CMI.h>

using namespace entropy::sparse;

double entropy::sparse::MC_W(ULContainer* W2, ULContainer* W1, ULContainer* A1, int mode)
{
  return CMI(W2, W1, A1, mode);
}
