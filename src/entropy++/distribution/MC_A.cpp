#include "MC_A.h"

#include <entropy++/distribution/CMI.h>

double entropy::distribution::MC_A(double*** pw2a1w1, int dimX, int dimY, int dimZ)
{
  return entropy::distribution::CMI(pw2a1w1, dimX, dimY, dimZ);
}
