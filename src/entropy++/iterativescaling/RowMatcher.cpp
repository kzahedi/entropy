#include "RowMatcher.h"

using namespace std;
using namespace entropy::iterativescaling;

RowMatcher::RowMatcher(int xSize)
{
  _x_rows = new Ivector[xSize];
  _y_rows = new Ivector[xSize];
  _rows   = xSize;

  for(int i = 0; i < xSize; i++) _x_rows[i].resize(0);
  for(int i = 0; i < xSize; i++) _y_rows[i].resize(0);
}

void RowMatcher::add_x(int delta_index, int x_index)
{
  _x_rows[delta_index].push_back(x_index);
}

void RowMatcher::add_y(int delta_index, int y_index)
{
  _y_rows[delta_index].push_back(y_index);
}

vector<int>::iterator RowMatcher::x_begin(int index)
{
  return _x_rows[index].begin();
}

vector<int>::iterator RowMatcher::x_end(int index)
{
  return _x_rows[index].end();
}

vector<int>::iterator RowMatcher::y_begin(int index)
{
  return _y_rows[index].begin();
}

vector<int>::iterator RowMatcher::y_end(int index)
{
  return _y_rows[index].end();
}
