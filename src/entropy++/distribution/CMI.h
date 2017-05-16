#ifndef __CMI_ON_DISTRIBUTION_H__
#define __CMI_ON_DISTRIBUTION_H__

namespace entropy
{
  namespace distribution
  {
    // MI(X;Y|Z)
    double CMI(double*** pxyz, int dimX, int dimY, int dimZ);
  }
}

#endif // __CMI_H__

