#include <entropy++/iterativescaling/Delta.h>

using namespace std;
using namespace entropy;
using namespace entropy::iterativescaling;

Delta::Delta(vector<unsigned long> xValues, vector<int> xIndices, vector<unsigned long> yValues, vector<int> yIndices)
{
  _xIndices = xIndices;
  _yIndices = yIndices;

  _xValues.resize(0);
  _yValues.resize(0);

  for(int i = 0; i < _xIndices.size(); i++)
  {
    _xValues.push_back(xValues[_xIndices[i]]);
  }

  for(int i = 0; i < _yIndices.size(); i++)
  {
    _yValues.push_back(yValues[_yIndices[i]]);
  }

  _observed = 0.0;
  _expected = 0.0;
  _lambda   = 1.0;
  _conditionalProbability = -1.0;
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

bool Delta::matchX(vector<unsigned long>& xValues)
{
  for(int i = 0; i < (int)_xIndices.size(); i++)
  {
    if(_xValues[i] != xValues[_xIndices[i]]) return false;
  }
  return true;
}

bool Delta::matchXY(vector<unsigned long>& xValues, vector<unsigned long>& yValues)
{
  for(int i = 0; i < (int)_xValues.size(); i++)
  {
    if(_xValues[i] != xValues[_xIndices[i]]) return false;
  }
  for(int i = 0; i < (int)_yValues.size(); i++)
  {
    if(_yValues[i] != yValues[_yIndices[i]]) return false;
  }
  return true;
}

bool Delta::matchP(vector<unsigned long>& xValues, vector<unsigned long>& yValues)
{
  if(xValues.size() != _xValues.size() ||
     yValues.size() != _yValues.size())
  {
    return false;
  }

  for(int i = 0; i < (int)xValues.size(); i++)
  {
    if(xValues[i] != _xValues[i]) return false;
  }

  for(int i = 0; i < (int)yValues.size(); i++)
  {
    if(yValues[i] != _yValues[i]) return false;
  }

  return false;
}
