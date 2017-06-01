#include <entropy++/iterativescaling/Delta.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace entropy;
using namespace entropy::iterativescaling;

Delta::Delta(vector<unsigned long> xv, vector<int> xi, vector<unsigned long> yv, vector<int> yi)
{

  _xIndices.resize(0);
  _yIndices.resize(0);
  _xValues.resize(0);
  _yValues.resize(0);

  for(int i = 0; i < xi.size();        i++) _xIndices.push_back(xi[i]);
  for(int i = 0; i < yi.size();        i++) _yIndices.push_back(yi[i]);
  for(int i = 0; i < _xIndices.size(); i++) _xValues.push_back(xv[_xIndices[i]]);
  for(int i = 0; i < _yIndices.size(); i++) _yValues.push_back(yv[_yIndices[i]]);

  _observed               = 1.0;
  _expected               = 0.0;
  _lambda                 = 1.0;
  _conditionalProbability = -1.0;
  _inputOnly              = false;
  _outputOnly             = false;
}

void Delta::incObserved()
{
  _observed += 1.0;
}

void Delta::setObserved(double v)
{
#ifdef USE_OPENMP
  _mutex.lock();
#endif // USE_OPENMP
  _observed = v;
#ifdef USE_OPENMP
  _mutex.lock();
#endif // USE_OPENMP
}

double Delta::observed()
{
  return _observed;
}

double Delta::expected()
{
  return _expected;
}

void Delta::updateExpected(double v)
{
#ifdef USE_OPENMP
  _mutex.lock();
#endif // USE_OPENMP
  _expected += v;
#ifdef USE_OPENMP
  _mutex.unlock();
#endif // USE_OPENMP
}

void Delta::setExpected(double v)
{
#ifdef USE_OPENMP
  _mutex.lock();
#endif // USE_OPENMP
  _expected = v;
#ifdef USE_OPENMP
  _mutex.unlock();
#endif // USE_OPENMP
}

void Delta::setLambda(double v)
{
#ifdef USE_OPENMP
  _mutex.lock();
#endif // USE_OPENMP
  _lambda = v;
#ifdef USE_OPENMP
  _mutex.unlock();
#endif // USE_OPENMP
}

void Delta::updateLambda(double v)
{
#ifdef USE_OPENMP
  _mutex.lock();
#endif // USE_OPENMP
  _lambda += v;
#ifdef USE_OPENMP
  _mutex.unlock();
#endif // USE_OPENMP
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

bool Delta::matchX(vector<unsigned long> xv)
{
  if(_outputOnly == true) return true;
  for(int i = 0; i < (int)_xValues.size(); i++) if(_xValues[i] != xv[_xIndices[i]]) return false;
  return true;
}

bool Delta::matchY(vector<unsigned long> yv)
{
  if(_inputOnly == true) return true;
  for(int i = 0; i < (int)_yValues.size(); i++) if(_yValues[i] != yv[_yIndices[i]]) return false;
  return true;
}

bool Delta::matchXY(vector<unsigned long> xv, vector<unsigned long> yv)
{
  if(_outputOnly == false)
    for(int i = 0; i < (int)_xValues.size(); i++) if(_xValues[i] != xv[_xIndices[i]]) return false;
  if(_inputOnly == false)
    for(int i = 0; i < (int)_yValues.size(); i++) if(_yValues[i] != yv[_yIndices[i]]) return false;
  return true;
}

bool Delta::match(vector<unsigned long> xv, vector<unsigned long> yv)
{
  if(_outputOnly == false && _inputOnly == false)
    if(xv.size() != _xValues.size() ||
       yv.size() != _yValues.size()) return false;
  if(_outputOnly == false && _inputOnly == true)
    if(xv.size() != _xValues.size()) return false;
  if(_outputOnly == true  && _inputOnly == false)
    if(yv.size() != _yValues.size()) return false;

  if(_outputOnly == false)
    for(int i = 0; i < (int)xv.size(); i++) if(xv[i] != _xValues[i]) return false;
  if(_inputOnly == false)
    for(int i = 0; i < (int)yv.size(); i++) if(yv[i] != _yValues[i]) return false;

  return true;
}

void Delta::setInputOnly()
{
  _inputOnly = true;
}

bool Delta::isInputOnly()
{
  return _inputOnly;
}

void Delta::setOutputOnly()
{
  _outputOnly = true;
}

bool Delta::isOutputOnly()
{
  return _outputOnly;
}
