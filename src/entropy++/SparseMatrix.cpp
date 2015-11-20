#include <entropy++/SparseMatrix.h>

#include <iostream>
#include <sstream>
#include <math.h>
#include <assert.h>

using namespace std;

SparseMatrix::SparseMatrix()
{
  _default   = 0.0;
  _size      = 0;
}

SparseMatrix::SparseMatrix(double def)
{
  _size      = 0;
  _default   = def;
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

void SparseMatrix::__set(const int f, const int s, int t, const double value)
{
  for(int i = 0; i < _indices.size(); i++)
  {
    MatrixIndex mi = _indices[i];
    if(mi.first == f && mi.second == s && mi.third == t)
    {
      _values[i] = value;
      return;
    }
  }
  MatrixIndex mi(f, s, t);
  _indices.push_back(mi);
  _values.push_back(value);
  if(f > _size) _size = f;
  if(s > _size) _size = s;
  if(t > _size) _size = t;
}

double SparseMatrix::__get(const int f, const int s, const int t) const
{
  double value = _default;
  for(int i = 0; i < _indices.size(); i++)
  {
    MatrixIndex mi = _indices[i];
    if(mi.first == f && mi.second == s && mi.third == t)
    {
      value = _values[i];
      break;
    }
  }
  return value;
}

double SparseMatrix::operator()(int first) const throw(MatrixException)
{
  return __get(first, -1, -1);
}

double& SparseMatrix::operator()(int first) throw(MatrixException)
{
  for(int i = 0; i < _indices.size(); i++)
  {
    MatrixIndex mi = _indices[i];
    if(mi.first == first)
    {
      return _values[i];
    }
  }
  MatrixIndex mi(first, -1, -1);
  _indices.push_back(mi);
  _values.push_back(_default);
  if(first  > _size) _size = first;
  return _values[_values.size() - 1];
}

double SparseMatrix::operator()(int first, int second) const throw(MatrixException)
{
  return __get(first, second, -1);
}

double& SparseMatrix::operator()(int first, int second) throw(MatrixException)
{
  for(int i = 0; i < _indices.size(); i++)
  {
    MatrixIndex mi = _indices[i];
    if(mi.first == first && mi.second == second)
    {
      return _values[i];
    }
  }
  MatrixIndex mi(first, second, -1);
  _indices.push_back(mi);
  _values.push_back(_default);
  if(first  > _size) _size = first;
  if(second > _size) _size = second;
  return _values[_values.size() - 1];
}

double SparseMatrix::operator()(int first, int second, int third) const throw(MatrixException)
{
  return __get(first, second, third);
}

double& SparseMatrix::operator()(int first, int second, int third) throw(MatrixException)
{
  for(int i = 0; i < _indices.size(); i++)
  {
    MatrixIndex mi = _indices[i];
    if(mi.first == first && mi.second == second)
    {
      return _values[i];
    }
  }
  MatrixIndex mi(first, second, third);
  _indices.push_back(mi);
  _values.push_back(_default);
  if(first  > _size) _size = first;
  if(second > _size) _size = second;
  if(third  > _size) _size = third;
  return _values[_values.size() - 1];
}

void SparseMatrix::__copy(const SparseMatrix &m)
{
  _indices.clear();
  _values.clear();
  for(int i = 0; i < (int)m._indices.size(); i++)
  {
    MatrixIndex mi = m._indices[i];
    double      v  = m._values[i];

    MatrixIndex nmi(mi.first, mi.second, mi.third);
    _indices.push_back(nmi);
    _values.push_back(v);
  }
  _size      = m._size;
}

SparseMatrix& SparseMatrix::operator+=(const SparseMatrix &m) throw(MatrixException)
{
  for(int i = 0; i < m._indices.size(); i++)
  {
    MatrixIndex mi       = m._indices[i];
    double value         = m._values[i];
    double current_value = __get(mi.first, mi.second, -1);
    __set(mi.first, mi.second, mi.third, current_value + value);
  }
  return *this;
}

SparseMatrix& SparseMatrix::operator-=(const SparseMatrix &m) throw(MatrixException)
{
  for(int i = 0; i < m._indices.size(); i++)
  {
    MatrixIndex mi       = m._indices[i];
    double value         = m._values[i];
    double current_value = __get(mi.first, mi.second, -1);
    __set(mi.first, mi.second, -1, current_value - value);
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
  for(int i = 0; i < _indices.size(); i++)
  {
    MatrixIndex mi       = _indices[i];
    double current_value = __get(mi.first, mi.second, -1);
    __set(mi.first, mi.second, -1, current_value - d);
  }
  return *this;
}

SparseMatrix& SparseMatrix::operator/=(const double d)
{
  for(int i = 0; i < _indices.size(); i++)
  {
    MatrixIndex mi       = _indices[i];
    double current_value = __get(mi.first, mi.second, -1);
    __set(mi.first, mi.second, -1, current_value / d);
  }
  return *this;
}

SparseMatrix operator*(double factor, const SparseMatrix& m) 
{ 
  SparseMatrix R;

  for(int i = 0; i < m.size(); i++)
  {
    MatrixIndex mi = m.getmi(i);
    int r = mi.first;
    int c = mi.second;
    cout << "row: " << r << " col: " << c << " value = " << m(r,c) << " -> ";
    R(r, c) = m(r,c) * factor;
    cout << R(r,c) << endl;
  }
  return R;
}

int SparseMatrix::size() const
{
  return _indices.size();
}

MatrixIndex SparseMatrix::getmi(int i) const
{
  return _indices[i];
}

double SparseMatrix::get(int index)
{
  return _values[index];
}
