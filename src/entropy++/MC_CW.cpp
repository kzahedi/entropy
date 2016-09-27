#include <entropy++/MC_CW.h>

#include <entropy++/MI.h>
#include <entropy++/ConditionalEntropy.h>

double entropy::MC_CW(ULContainer* W2,
                      ULContainer* W1,
                      ULContainer* A1,
                      int mc_cw_mode,
                      int mode) throw (EntropyException)
{
  switch(mc_cw_mode)
  {
    case MC_CW_MODE_ENTROPY:
      return 0.0;
      return entropy::ConditionalEntropy(W2, A1)
        - entropy::ConditionalEntropy(W2, W1);
      break;
    case MC_CW_MODE_MI:
      return entropy::MI(W2, W1, mode)
        - entropy::MI(W2, A1, mode);
      break;
    default:
      throw EntropyException("Unknown MC_CW mode");
  }

}
