#include <entropy++/sparse/MC_CW.h>
#include <entropy++/sparse/MI.h>

using namespace entropy::sparse;

double entropy::sparse::MC_CW(ULContainer* W2, ULContainer* W1, ULContainer* A1, int mode)
{
  double mi_w2_w1 = entropy::sparse::MI(W2, W1, mode);
  double mi_w2_a1 = entropy::sparse::MI(W2, A1, mode);
  return mi_w2_w1 - mi_w2_a1;
}
