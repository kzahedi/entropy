#ifndef __CMI_H__
#define __CMI_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

namespace entropy
{
  /*! \brief Conditional Mutual Information
   *
   *  This function returns I(X;Y|Z)
   *  \param mode currently only supports EMPERICAL, which is a histogram based
   *  method to estimate the probability distributions
   */
  double CMI(ULContainer* X, ULContainer* Y, ULContainer* Z, int mode = EMPERICAL);
}

#endif // __CMI_H__
