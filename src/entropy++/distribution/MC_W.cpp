#include "MC_W.h"

#include <entropy++/distribution/CMI.h>

double entropy::distribution::MC_W(double*** pw2w1a1, int dimX, int dimY, int dimZ)
{
  return entropy::distribution::CMI(pw2w1a1, dimX, dimY, dimZ);
}
