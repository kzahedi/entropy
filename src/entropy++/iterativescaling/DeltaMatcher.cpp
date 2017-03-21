#include "DeltaMatcher.h"
#include "Delta.h"

using namespace std;
using namespace entropy::iterativescaling;

DeltaMatcher::DeltaMatcher(int xSize)
{
  _deltas = new Dvector[xSize];
  _rows   = xSize;

  for(int i = 0; i < xSize; i++) _deltas[i].resize(0);
}

void DeltaMatcher::add(int index, Delta* d)
{
  _deltas[index].push_back(d);
}

vector<Delta*>::iterator DeltaMatcher::begin(int index)
{
  return _deltas[index].begin();
}

vector<Delta*>::iterator DeltaMatcher::end(int index)
{
  return _deltas[index].end();
}
