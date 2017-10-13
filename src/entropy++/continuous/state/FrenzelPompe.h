#ifndef __FRENZEL_POMPE_H__
#define __FRENZEL_POMPE_H__

#include <entropy++/continuous/Functions.h>
#include <entropy++/Container.h>
#include <entropy++/EntropyException.h>

namespace entropy
{
  namespace continuous
  {
    namespace state
    {
      // Calculates the conditional mutual information on continuous data,
      // according to S. Frenzel and B. Pompe.
      // Partial mutual information for coupling analysis
      // of multivariate time series.
      // Phys. Rev. Lett., 99:204101, Nov 2007.
      double FrenzelPompe(entropy::DContainer* xyz,
                          vector<int> xIndices,
                          vector<int> yIndices,
                          vector<int> zIndices,
                          int k);
    }
  }
}

#endif // __FRENZEL_POMPE_H__
