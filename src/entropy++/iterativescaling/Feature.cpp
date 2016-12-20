#include <entropy++/iterativescaling/Feature.h>

using namespace entropy;
using namespace entropy::iterativescaling;


Feature::Feature(int xListIndex, int yListIndex)
{
  _xListIndex = xListIndex;
  _yListIndex = yListIndex;
}

int Feature::xListIndex()
{
  return _xListIndex;
}

int Feature::yListIndex()
{
  return _yListIndex;
}
