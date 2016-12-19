#include <entropy++/iterativescaling/Feature.h>

using namespace entropy;
using namespace entropy::iterativescaling;


Feature::Feature(int xListIndex, int yListIndex)
{
  _xListIndex             = xListIndex;
  _yListIndex             = yListIndex;
  _xAlphabetSize          = 0.0;
  _yAlphabetSize          = 0.0;
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

void Feature::setXAlphabetSize(double x)
{
  _xAlphabetSize = x;
}

double Feature::getXAlphabetSize()
{
  return _xAlphabetSize;
}

void Feature::setYAlphabetSize(double y)
{
  _yAlphabetSize = y;
}

double Feature::getYAlphabetSize()
{
  return _yAlphabetSize;
}


