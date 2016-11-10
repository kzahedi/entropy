#include <entropy++/sparse/MC_CW.h>
#include <entropy++/sparse/MI.h>
#include <entropy++/sparse/ConditionalEntropy.h>

using namespace entropy::sparse;

double entropy::sparse::MC_CW(ULContainer* W2,
                              ULContainer* W1,
                              ULContainer* A1,
                              int mc_cw_mode,
                              int mode)
  throw (EntropyException)
{
  switch(mc_cw_mode)
  {
    case MC_CW_MODE_ENTROPY:
      return entropy::sparse::ConditionalEntropy(W2, A1)
        - entropy::sparse::ConditionalEntropy(W2, W1);
      break;
    case MC_CW_MODE_MI:
      return entropy::sparse::MI(W2, W1, mode)
        - entropy::sparse::MI(W2, A1, mode);
      break;
    default:
      throw EntropyException("Unknown MC_CW mode");
  }
}
