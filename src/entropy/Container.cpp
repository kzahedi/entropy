#include <entropy/Container.h>

#include <iostream>
#include <assert.h>

# define MIN(a,b) (((a)<(b))?a:b)

using namespace std;

Container::Container(int rows, int columns)
{
  _data         = new double*[rows];
  _domains      = new double*[columns];
  _bins         = new int[columns];
  _mode         = CONTAINER_DISCRETISE_UNIFORM;

  _domainsGiven = false;
  _binsGiven    = false;

  for(int i = 0; i < columns; i++)
  {
    _domains[i] = new double[2];
  }

  for(int i = 0; i < rows; i++)
  {
    _data[i] = new double[columns];
  }

  for(int r = 0; r < rows; r++)
  {
    for(int c = 0; c < columns; c++)
    {
      _data[r][c] = 0.0;
    }
  }

  _rows      = rows;
  _columns   = columns;
  _fillIndex = 0;
}

Container::~Container()
{
  for(int r = 0; r < _rows;    r++) delete _data[r];
  for(int c = 0; c < _columns; c++) delete _domains[c];

  delete _data;
  delete _domains;
  delete _bins;
}

double Container::operator()(const int row, const int column) const
{
  assert(row    < _rows);
  assert(column < _columns);
  return _data[row][column];
}

const Container& Container::operator<<(const double &value) const
{
  assert(_fillIndex < _rows * _columns);

  int c = _fillIndex % _columns;
  int r = _fillIndex / _columns;

  _data[r][c] = value;
  Container *co = (Container*)this;
  co->_fillIndex++;

  return *this;
}

int Container::rows()
{
  return _rows;
}

int Container::columns()
{
  return _columns;
}

double Container::get(int row, int column)
{
  assert(row    < _rows);
  assert(column < _columns);
  return _data[row][column];
}

void Container::setBinSizes(int *bins)
{
  for(int c = 0; c < _columns; c++)
  {
    _bins[c] = bins[c];
  }
  _binsGiven = true;
}

void Container::setDomains(double **domains)
{
  for(int c = 0; c < _columns; c++)
  {
    _domains[c][0] = domains[c][0];
    _domains[c][1] = domains[c][1];
  }
  _domainsGiven = true;
}


Container* Container::discretise()
{
  switch(_mode)
  {
    case CONTAINER_DISCRETISE_UNIFORM:
      return __uniformDiscretisation();
      break;
    default:
      cerr << "Unknown discretisation mode given: " << _mode << endl;
      break;
  }
  return NULL;
}

void Container::setDiscretisationMode(int mode)
{
  _mode = mode;
}

int Container::__discretiseAndCombineValues(double *values)
{
  assert(_domainsGiven && _binsGiven);

  int value  = 0;
  int factor = 1;

  for(int c = 0; c < _columns; c++)
  {
    double mapped = ((values[c]      - _domains[c][0]) /
                     (_domains[c][1] - _domains[c][0]) * _bins[c]);
    double cropped = MIN(_bins[c]-1, mapped);
    if(c > 0) factor *= _bins[c-1];
    value += (int)(cropped * factor);
  }

  return value;
}

Container* Container::__uniformDiscretisation()
{
  Container *copy = new Container(_rows, 1);
  for(int r = 0; r < _rows; r++)
  {
    double v = (double)__discretiseAndCombineValues(_data[r]);
    (*copy) << v;
  }
  return copy;
}
