#include <entropy++/iterativescaling/Feature.h>

using namespace entropy;
using namespace entropy::iterativescaling;


Feature::Feature(int xListIndex, int yListIndex)
{
  _xListIndex = xListIndex;
  _yListIndex = yListIndex;
  _alphabetSize = 0.0;
}

int Feature::xListIndex()
{
  return _xListIndex;
}

int Feature::yListIndex()
{
  return _yListIndex;
}

// obs(x)
void Feature::setUniqueXCount(int index, int count)
{
  if(_uniqueXCount.size() < index + 1)
  {
    _uniqueXCount.resize(index+1);
  }
  _uniqueXCount[index] = count;
}


int Feature::getUniqueXCount(int index)
{
  return _uniqueXCount[index];
}

void Feature::setRemainingAlphabetSize(double s)
{
  _alphabetSize = s;
}

double Feature::getRemainingAlphabetSize()
{
  return _alphabetSize;
}

void Feature::setYAlphabetSize(double y)
{
  _yAlphabetSize = y;
}

double Feature::yAlphabetSize()
{
  return _yAlphabetSize;
}
