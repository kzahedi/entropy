#include <entropy++/SparseMatrix.h>

#include <iostream>
#include <sstream>
#include <math.h>

using namespace std;

SparseMatrix::SparseMatrix()
{
  _rows = 0;
  _cols = 0;
}

SparseMatrix::SparseMatrix(SparseMatrix &m)
{
  __copy(m);
}

SparseMatrix::SparseMatrix(const SparseMatrix &m)
{
  __copy(m);
}

SparseMatrix::~SparseMatrix()
{
}

void SparseMatrix::reset()
{
  _indices.clear();
  _values.clear();
}

void SparseMatrix::__set(const int row, const int col, const double value)
{
  for(int i = 0; i < _indices.size(); i++)
  {
    matrixindex mi = _indices[i];
    if(mi.first == row && mi.second == col)
    {
      _values[i] = value;
      return;
    }
  }
  matrixindex mi = make_pair(row, col);
  _indices.push_back(mi);
  _values.push_back(value);
}

double SparseMatrix::__get(const int row, const int col) const
{
  double value = 0.0;
  for(int i = 0; i < _indices.size(); i++)
  {
    matrixindex mi = _indices[i];
    if(mi.first == row && mi.second == col)
    {
      value = _values[i];
      break;
    }
  }
  return value;
}

double SparseMatrix::operator()(int row, int col) const throw(MatrixException)
{
  return __get(row, col);
}

double& SparseMatrix::operator()(int row, int col) throw(MatrixException)
{
  for(int i = 0; i < _indices.size(); i++)
  {
    matrixindex mi = _indices[i];
    if(mi.first == row && mi.second == col)
    {
      return _values[i];
    }
  }
  matrixindex mi = make_pair(row, col);
  _indices.push_back(mi);
  _values.push_back(0.0);
  return _values[_values.size() - 1];
}

double SparseMatrix::get(int row, int col)
{
  return __get(row, col);
}

void SparseMatrix::__copy(const SparseMatrix &m)
{
  for(int i = 0; i < m._indices.size(); i++)
  {
    matrixindex mi = m._indices[i];
    double      v  = m._values[i];

    matrixindex nmi = make_pair(mi.first, mi.second);
    _indices.push_back(nmi);
    _values.push_back(v);
  }
}

SparseMatrix& SparseMatrix::operator+=(const SparseMatrix &m) throw(MatrixException)
{
  for(int i = 0; m._indices.size(); i++)
  {
    matrixindex mi       = m._indices[i];
    double value         = m._values[i];
    double current_value = __get(mi.first, mi.second);
    __set(mi.first, mi.second, current_value + value);
  }
  return *this;
}

SparseMatrix& SparseMatrix::operator-=(const SparseMatrix &m) throw(MatrixException)
{
  for(int i = 0; m._indices.size(); i++)
  {
    matrixindex mi       = m._indices[i];
    double value         = m._values[i];
    double current_value = __get(mi.first, mi.second);
    __set(mi.first, mi.second, current_value - value);
  }
  return *this;
}

SparseMatrix& SparseMatrix::operator*=(const double factor)
{
  for(int i = 0; i < _values.size(); i++)
  {
    _values[i] *= factor;
  }
  return *this;
}

const SparseMatrix SparseMatrix::operator+(const SparseMatrix &m)
{
  SparseMatrix r = *this;
  r += m;
  return r;
}

const SparseMatrix SparseMatrix::operator-(const SparseMatrix &m)
{
  SparseMatrix r = *this;
  r -= m;
  return r;
}

const SparseMatrix SparseMatrix::operator*(const double factor)
{
  SparseMatrix r = *this;
  r *= factor;
  return r;
}

SparseMatrix& SparseMatrix::operator=(const SparseMatrix &m)
{
  if(this == &m) return *this;
  __copy(m);
  return *this;
}

SparseMatrix& SparseMatrix::operator-=(const double d)
{
  for(int i = 0; _indices.size(); i++)
  {
    matrixindex mi       = _indices[i];
    double current_value = __get(mi.first, mi.second);
    __set(mi.first, mi.second, current_value - d);
  }
  return *this;
}

SparseMatrix& SparseMatrix::operator/=(const double d)
{
  for(int i = 0; _indices.size(); i++)
  {
    matrixindex mi       = _indices[i];
    double current_value = __get(mi.first, mi.second);
    __set(mi.first, mi.second, current_value / d);
  }
  return *this;
}

SparseMatrix& SparseMatrix::operator*=(const SparseMatrix &m)
{
  for(int i = 0; m._indices.size(); i++)
  {
    matrixindex mi       = m._indices[i];
    double value         = m._values[i];
    double current_value = __get(mi.first, mi.second);
    __set(mi.first, mi.second, current_value * value);
  }
  return *this;
}

SparseMatrix operator*(double factor, const SparseMatrix& m) 
{ 
  SparseMatrix R;

  for(int i = 0; i < m.size(); i++)
  {
    matrixindex mi = m.getmi(i);
    int r = mi.first;
    int c = mi.second;
    R(r, c) = m(r,c) * factor;
  }
  return R;
}

int SparseMatrix::size() const
{
  return _indices.size();
}

matrixindex SparseMatrix::getmi(int i) const
{
  return _indices[i];
}


