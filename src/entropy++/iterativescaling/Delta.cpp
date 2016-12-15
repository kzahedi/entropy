#include <entropy++/iterativescaling/Delta.h>

using namespace std;
using namespace entropy;
using namespace entropy::iterativescaling;

Delta::Delta(int xUniqueIndex, int yUniqueIndex)
{
  _xUniqueIndex           = xUniqueIndex;
  _yUniqueIndex           = yUniqueIndex;
  _observed               = 0.0;
  _expected               = 0.0;
  _lambda                 = 1.0;
  _conditionalProbability = -1.0;
}

bool Delta::match(int xUniqueIndex, int yUniqueIndex)
{
  return (_xUniqueIndex == xUniqueIndex &&
          _yUniqueIndex == yUniqueIndex);
}

int Delta::getUniqueXIndex()
{
  return _xUniqueIndex;
}

int Delta::getUniqueYIndex()
{
  return _yUniqueIndex;
}

void Delta::incObserved()
{
  _observed += 1.0;
}

void Delta::setObserved(double v)
{
  _observed = v;
}

double Delta::observed()
{
  return _observed;
}

double Delta::expected()
{
  return _expected;
}

void Delta::setExpected(double v)
{
  _expected = v;
}

void Delta::setLambda(double v)
{
  _lambda = v;
}

double Delta::lambda()
{
  return _lambda;
}

void Delta::setConditionalProbability(double value)
{
  _conditionalProbability = value;
}

double Delta::conditionalProbability()
{
  return _conditionalProbability;
}

void Delta::setMarginalProbability(double value)
{
  _marginalProbability = value;
}

double Delta::marginalProbability()
{
  return _marginalProbability;
}

