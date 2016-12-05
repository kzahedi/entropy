#include <entropy++/SparseMatrix.h>

#include <iostream>
#include <sstream>
#include <math.h>
#include <assert.h>

using namespace std;

entropy::SparseMatrix::SparseMatrix()
{
  _default = 0.0;
  _size    = 0;
}

entropy::SparseMatrix::SparseMatrix(double def)
{
  _size    = 0;
  _default = def;
}

entropy::SparseMatrix::SparseMatrix(SparseMatrix &m)
{
  __copy(m);
}

entropy::SparseMatrix::SparseMatrix(const SparseMatrix &m)
{
  __copy(m);
}

entropy::SparseMatrix::~SparseMatrix()
{ }

void entropy::SparseMatrix::reset()
{
  _indices.clear();
  _values.clear();
}

void entropy::SparseMatrix::__set(const int f, const int s, const int t, const double value)
{
  for(int i = 0; i < (int)_indices.size(); i++)
  {
    entropy::MatrixIndex mi = _indices[i];
    if(mi.first == f && mi.second == s && mi.third == t)
    {
      _values[i] = value;
      return;
    }
  }

  entropy::MatrixIndex mi(f, s, t);
  _indices.push_back(mi);
  _values.push_back(value);
  if(f > _size) _size = f;
  if(s > _size) _size = s;
  if(t > _size) _size = t;
}

double entropy::SparseMatrix::__get(const int f, const int s, const int t) const
{
  for(int i = 0; i < (int)_indices.size(); i++)
  {
    entropy::MatrixIndex mi = _indices[i];
    if(mi.first == f && mi.second == s && mi.third == t)
    {
      // if(mi.first == 0 && mi.second == 0 && mi.third == 1)
      // {
      // }
      // if(f == 0 && s == 0 && t == 1)
      // {
      // }
      return _values[i];
    }
  }
  // if(f == 0 && s == 0 && t == 1)
  // {
  // }
  return _default;
}

double entropy::SparseMatrix::operator()(int first) const
  throw(EntropyException)
{
  return __get(first, -1, -1);
}

double& entropy::SparseMatrix::operator()(int first)
  throw(EntropyException)
{
  for(int i = 0; i < (int)_indices.size(); i++)
  {
    entropy::MatrixIndex mi = _indices[i];
    if(mi.first == first)
    {
      return _values[i];
    }
  }
  entropy::MatrixIndex mi(first, -1, -1);
  _indices.push_back(mi);
  _values.push_back(_default);
  if(first  > _size) _size = first;
  return _values[_values.size() - 1];
}

double entropy::SparseMatrix::operator()(int first, int second) const
  throw(EntropyException)
{
  return __get(first, second, -1);
}

double& entropy::SparseMatrix::operator()(int first, int second)
  throw(EntropyException)
{
  for(int i = 0; i < (int)_indices.size(); i++)
  {
    entropy::MatrixIndex mi = _indices[i];
    if(mi.first == first && mi.second == second)
    {
      return _values[i];
    }
  }
  entropy::MatrixIndex mi(first, second, -1);
  _indices.push_back(mi);
  _values.push_back(_default);
  if(first  > _size) _size = first;
  if(second > _size) _size = second;
  return _values[_values.size() - 1];
}

double entropy::SparseMatrix::operator()(int first, int second, int third) const
  throw(EntropyException)
{
  return __get(first, second, third);
}

double& entropy::SparseMatrix::operator()(int first, int second, int third)
  throw(EntropyException)
{
  for(int i = 0; i < (int)_indices.size(); i++)
  {
    entropy::MatrixIndex mi = _indices[i];
    if(mi.first == first && mi.second == second && mi.third == third)
    {
      return _values[i];
    }
  }
  entropy::MatrixIndex mi(first, second, third);
  _indices.push_back(mi);
  _values.push_back(_default);
  if(first  > _size) _size = first;
  if(second > _size) _size = second;
  if(third  > _size) _size = third;
  return _values[_values.size() - 1];
}

void entropy::SparseMatrix::__copy(const SparseMatrix &m)
{
  _indices.clear();
  _values.clear();
  for(int i = 0; i < (int)m._indices.size(); i++)
  {
    entropy::MatrixIndex mi = m._indices[i];
    double      v  = m._values[i];

    entropy::MatrixIndex nmi(mi.first, mi.second, mi.third);
    _indices.push_back(nmi);
    _values.push_back(v);
  }
  _size    = m._size;
  _default = m._default;
}

entropy::SparseMatrix& entropy::SparseMatrix::operator+=(const SparseMatrix &m)
  throw(EntropyException)
{
  for(int i = 0; i < (int)m._indices.size(); i++)
  {
    entropy::MatrixIndex mi       = m._indices[i];
    double value         = m._values[i];
    double current_value = __get(mi.first, mi.second, mi.third);
    __set(mi.first, mi.second, mi.third, current_value + value);
  }
  return *this;
}

entropy::SparseMatrix& entropy::SparseMatrix::operator-=(const SparseMatrix &m)
  throw(EntropyException)
{
  for(int i = 0; i < (int)m._indices.size(); i++)
  {
    entropy::MatrixIndex mi       = m._indices[i];
    double value         = m._values[i];
    double current_value = __get(mi.first, mi.second, mi.third);
    __set(mi.first, mi.second, mi.third, current_value - value);
  }
  return *this;
}

entropy::SparseMatrix& entropy::SparseMatrix::operator*=(const double factor)
{
  for(int i = 0; i < (int)_values.size(); i++)
  {
    _values[i] *= factor;
  }
  return *this;
}

const entropy::SparseMatrix entropy::SparseMatrix::operator+(const SparseMatrix &m)
{
  SparseMatrix r = *this;
  r += m;
  return r;
}

const entropy::SparseMatrix entropy::SparseMatrix::operator-(const SparseMatrix &m)
{
  SparseMatrix r = *this;
  r -= m;
  return r;
}

const entropy::SparseMatrix entropy::SparseMatrix::operator*(const double factor)
{
  SparseMatrix r = *this;
  r *= factor;
  return r;
}

entropy::SparseMatrix& entropy::SparseMatrix::operator=(const SparseMatrix &m)
{
  if(this == &m) return *this;
  __copy(m);
  return *this;
}

entropy::SparseMatrix& entropy::SparseMatrix::operator-=(const double d)
{
  for(int i = 0; i < (int)_indices.size(); i++)
  {
    entropy::MatrixIndex mi       = _indices[i];
    double current_value = __get(mi.first, mi.second, mi.third);
    __set(mi.first, mi.second, mi.third, current_value - d);
  }
  return *this;
}

entropy::SparseMatrix& entropy::SparseMatrix::operator/=(const double d)
{
  for(int i = 0; i < (int)_indices.size(); i++)
  {
    entropy::MatrixIndex mi       = _indices[i];
    double current_value = __get(mi.first, mi.second, mi.third);
    if(mi.first == 0 && mi.second == 0 && mi.third == 1)
    {
    }
    __set(mi.first, mi.second, mi.third, current_value / d);
  }
  return *this;
}

entropy::SparseMatrix operator*(double factor, const entropy::SparseMatrix& m) 
{ 
  entropy::SparseMatrix R;

  for(int i = 0; i < m.size(); i++)
  {
    entropy::MatrixIndex mi = m.getmi(i);
    int r = mi.first;
    int c = mi.second;
    R(r, c) = m(r,c) * factor;
  }
  return R;
}

int entropy::SparseMatrix::size() const
{
  return _indices.size();
}

entropy::MatrixIndex entropy::SparseMatrix::getmi(int i) const
{
  return _indices[i];
}

double entropy::SparseMatrix::get(int index)
{
  return _values[index];
}

double entropy::SparseMatrix::defaultValue()
{
  return _default;
}

bool entropy::SparseMatrix::available(int first, int second, int third)
{
  for(int i = 0; i < (int)_indices.size(); i++)
  {
    if(_indices[i].first  == first  &&
       _indices[i].second == second &&
       _indices[i].third  == third)
      return true;
  }
  return false;
}

bool entropy::SparseMatrix::available(int first, int second)
{
  return available(first, second, -1);
}

bool entropy::SparseMatrix::available(int first)
{
  return available(first, -1, -1);
}

double entropy::SparseMatrix::value(int first, int second, int third)
{
  return __get(first, second, third);
}

double entropy::SparseMatrix::value(int first, int second)
{
  return __get(first, second, -1);
}

double entropy::SparseMatrix::value(int first)
{
  return __get(first, -1, -1);
}

