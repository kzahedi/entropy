#include <entropy++/Container.h>

#include <iostream>
#include <assert.h>
#include <vector>

# define MIN(a,b) (((a)<(b))?a:b)

using namespace std;

Container::Container(int rows, int columns)
{
  _data         = new double*[rows];
  _domains      = new double*[columns];
  _bins         = new int[columns];
  _mode         = UNIFORM;

  _domainsGiven = false;
  _binsGiven    = false;
  _discretised  = false;

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

// Container(const Container);
Container& Container::operator=(const Container &c)
{
  if(_data != NULL)
  {
    for(int r = 0; r < _rows;    r++) delete _data[r];
    delete _data;
  }

  if(_domains != NULL)
  {
    for(int c = 0; c < _columns; c++) delete _domains[c];
    delete _domains;
  }

  if(_bins != NULL)
  {
    delete _bins;
  }

  this->_rows    = c._rows;
  this->_columns = c._columns;

  this->_data         = new double*[this->_rows];
  this->_domains      = new double*[this->_columns];
  this->_bins         = new int[this->_columns];
  this->_mode         = c._mode;

  this->_domainsGiven = c._domainsGiven;
  this->_binsGiven    = c._binsGiven;
  this->_discretised  = c._discretised;

  if(this->_domainsGiven)
  {
    for(int i = 0; i < this->_columns; i++)
    {
      this->_domains[i] = new double[2];
      for(int j = 0; j < 2; j++)
      {
        this->_domains[i][j] = c._domains[i][j];
      }
    }
  }

  if(this->_binsGiven)
  {
    for(int i = 0; i < this->_columns; i++)
    {
      this->_bins[i] = c._bins[i];
    }
  }

  for(int i = 0; i < this->_rows; i++)
  {
    this->_data[i] = new double[this->_columns];
    for(int j = 0; j < this->_columns; j++)
    {
      this->_data[i][j] = c._data[i][j];
    }
  }
  return *this;
}

Container::~Container()
{
  if(_data != NULL)
  {
    for(int r = 0; r < _rows;    r++) delete _data[r];
    delete _data;
  }

  if(_domains != NULL && _domainsGiven == true)
  {
    for(int c = 0; c < _columns; c++) delete _domains[c];
    delete _domains;
  }

  if(_bins != NULL && _binsGiven == true)
  {
    delete _bins;
  }

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

void Container::set(int row, int column, double value)
{
  assert(row    < _rows);
  assert(column < _columns);
  _data[row][column] = value;
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
  Container *dis = NULL;
  Container *com = NULL;
  switch(_mode)
  {
    case UNIFORM:
      dis = discretiseByColumn();
      com = dis->combineDiscretisedColumns();
      delete dis;
      return com;
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

void Container::__strip()
{
  vector<int> values;

  for(int r = 0; r < _rows; r++)
  {
    int  value = (int)_data[r][0];
    bool found = false;
    for(int i = 0; i < (int)values.size(); i++)
    {
      if(value == values[i])
      {
        found = true;
        break;
      }
    }
    if(found == false)
    {
      values.push_back(value);
    }
  }

  for(int r = 0; r < _rows; r++)
  {
    int value = (int)_data[r][0];
    for(int i = 0; i < (int)values.size(); i++)
    {
      if(value == values[i])
      {
        _data[r][0] = (double)i;
      }
    }
  }
}

Container* Container::__uniformDiscretisationByColumn()
{
  Container *copy = new Container(_rows, _columns);
  __copyProperties(copy);
  for(int c = 0; c < _columns; c++)
  {
    vector<int> values;
    for(int r = 0; r < _rows; r++)
    {
      assert(_domains[c][0] <= get(r,c) && get(r,c) <= _domains[c][1]);
      int mapped  = (int)(((get(r,c)       - _domains[c][0])
                         / (_domains[c][1] - _domains[c][0]))
                         * _bins[c]);
      int cropped = (int)MIN(_bins[c]-1, mapped);
      copy->set(r, c, cropped);
    }
  }
  copy->__strip();
  copy->_discretised = true;
  return copy;
}

bool Container::isDiscretised()
{
  return _discretised;
}

double& Container::operator()(const int row, const int column)
{
  assert(row    < _rows);
  assert(column < _columns);
  return _data[row][column];
}

Container* Container::drop(int n)
{
  if(n > 0) return __dropFirst(n);
  if(n < 0) return __dropLast(n);
  return copy();
}

Container* Container::__dropFirst(int n)
{
  Container *copy = new Container(_rows-abs(n), _columns);
  __copyProperties(copy);

  for(int i = 0; i < _rows-abs(n); i++)
  {
    for(int j = 0; j < _columns; j++)
    {
      (*copy) << _data[i+abs(n)][j];
    }
  }
  return copy;
}

Container* Container::__dropLast(int n)
{
  Container *copy = new Container(_rows-abs(n), _columns);
  __copyProperties(copy);

  for(int i = 0; i < _rows - abs(n); i++)
  {
    for(int j = 0; j < _columns; j++)
    {
      (*copy) << _data[i][j];
    }
  }
  return copy;
}

double Container::max()
{
  double m = max(0);
  for(int i = 1; i < _columns; i++)
  {
    double n = max(i);
    if(n > m) m = n;
  }
  return m;
}

double Container::max(int column)
{
  double m = _data[0][column];
  for(int i = 1; i < _rows; i++)
  {
    if(_data[i][column] > m) m = _data[i][column];
  }
  return m;
}

double Container::min()
{
  double m = min(0);
  for(int i = 1; i < _columns; i++)
  {
    double n = min(i);
    if(n < m) m = n;
  }
  return m;
}

double Container::min(int column)
{
  double m = _data[0][column];
  for(int i = 1; i < _rows; i++)
  {
    if(_data[i][column] < m) m = _data[i][column];
  }
  return m;
}


Container* Container::columns(int n, ...)
{
  vector<int> indices;
  va_list ap;
  va_start(ap, n);
  for(int i = 0; i < n; i++)
  {
     indices.push_back(va_arg(ap, int));
  }
  va_end(ap);

  Container *extracted = new Container(this->rows(), n);

  for(int r = 0; r < _rows; r++)
  {
    for(int i = 0; i < (int)indices.size(); i++)
    {
      (*extracted) << get(r, indices[i]);
    }
  }

  return extracted;
}

void Container::normaliseColumn(int c, double min, double max)
{
  for(int r = 0; r < _rows; r++)
  {
    _data[r][c] = (_data[r][c] - min) / (max - min);
  }
}

Container* Container::copy()
{
  Container *copy = new Container(_rows, _columns);
  __copyProperties(copy);

  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _columns; c++)
    {
      (*copy) << _data[r][c];
    }
  }
  return copy;
}

Container* Container::discretiseByColumn()
{
  switch(_mode)
  {
    case UNIFORM:
      return __uniformDiscretisationByColumn();
      break;
    default:
      cerr << "Unknown discretisation mode given: " << _mode << endl;
      break;
  }
  return NULL;
}

void Container::__copyProperties(Container* dst)
{
  dst->_binsGiven    = _binsGiven;
  dst->_domainsGiven = _domainsGiven;
  dst->_discretised  = _discretised;
  if(_binsGiven == true && dst->_columns == _columns)
  {
    for(int i = 0; i < _columns; i++)
    {
      dst->_bins[i] = _bins[i];
    }
  }
  if(_domainsGiven == true && dst->_columns == _columns)
  {
    for(int i = 0; i < _columns; i++)
    {
      dst->_domains[i][0] = _domains[i][0];
      dst->_domains[i][1] = _domains[i][1];
    }
  }
}

Container* Container::combineDiscretisedColumns()
{
  assert(_discretised == true);
  Container *copy = new Container(_rows, 1);
  __copyProperties(copy);

  for(int r = 0; r < _rows; r++)
  {
    int v = 0;
    int f = 1;
    for(int c = 0; c < _columns; c++)
    {
      if(c > 0) f = f * _bins[c-1];
      v += f * get(r,c);
    }
    *copy << v;
  }
  copy->__strip();
  return copy;
}

double Container::columnSum(int c)
{
  double s = 0.0;
  for(int r = 0; r < _rows; r++) s += _data[r][c];
  return s;
}
