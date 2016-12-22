#include <entropy++/Matrix.h>

#include <iostream>
#include <sstream>
#include <math.h>

using namespace std;
using namespace entropy;

Matrix::Matrix()
{
  _cell = NULL;
  _rows = 0;
  _cols = 0;
}

Matrix::Matrix(Matrix &m)
{
  __init(m.rows(), m.cols(), 0.0);
  __copy(m);
}

Matrix::Matrix(const Matrix &m)
{
  __init(m.rows(), m.cols(), 0.0);
  __copy(m);
}

Matrix::Matrix(int rows, int cols)
{
  _cell = NULL;
  __init(rows,cols,0.0); // standard matrix
}

Matrix::Matrix(int rows, int cols, double initialValue)
{
  __init(rows, cols, initialValue); // standard matrix
}

Matrix::Matrix(int rows, int cols, std::vector<double> initialValues)
{
  __init(rows,cols, 0.0);
  for(int i = 0; i < _rows; i++)
  {
    for(int j = 0; j < _cols; j++)
    {
      if((unsigned int)(i*_rows + j) > initialValues.size())
      {
        _cell[i][j] = 0.0;
      }
      else
      {
        _cell[i][j] = initialValues[i*_rows + j];
      }
    }
  }
}

Matrix::~Matrix()
{
  __deleteCells();
}

void Matrix::__init(const int rows, const int cols, double initialValue)
{
  if(rows < 0 || cols < 0)
  {
    std::stringstream oss;
    oss << "Invalid matrix initialisation: ";
    if(rows < 0)
    {
      oss << "The number of rows is negative: " << rows;
    }
    if(cols < 0)
    {
      if (rows < 0)
      {
        oss << " and the ";
      }
      else
      {
        oss << "The ";
      }
      oss << "number of columns is negative: " << cols;
    }
    throw EntropyException(oss.str());
  }

  _rows = rows;
  _cols = cols;

  _cell = new double* [rows];
  for(int i=0; i<rows; i++)
  {
    _cell[i] = new double [cols];
  }

  for(int i = 0; i < rows; i++)
  {
    for(int j = 0; j < cols; j++)
    {
      _cell[i][j] = initialValue;
    }
  }
}


double Matrix::operator()(int row, int col) const throw(EntropyException)
{
  if(row >= _rows || col >= _cols || row < 0 || col < 0)
  {
    std::stringstream oss;
    oss << "Invalid matrix index: ";
    if(row > _rows)
    {
      oss << "row too large: " << row << " > " << _rows;
    }
    if(row < 0)
    {
      oss << "row is negative: " << row;
    }
    if(col > _cols)
    {
      oss << "col too large: " << col << " > " << _cols;
    }
    if(col < 0)
    {
      oss << "col is negative: " << col;
    }
    throw EntropyException(oss.str());
  }
  return _cell[row][col];
}

double& Matrix::operator()(int row, int col) throw(EntropyException)
{
  __check(row,col);
  return _cell[row][col];
}

void Matrix::__set(const int row, const int col, const double value)
{
  __check(row,col);
  _cell[row][col] = value;
}

double Matrix::get(int row, int col)
{
  return _cell[row][col];
}

double Matrix::rowSum(int row)
{
  double sum = 0;
  for(int c = 0; c < _cols; c++)
  {
    sum += _cell[row][c];
  }
  return sum;
}

double Matrix::colSum(int col)
{
  double sum = 0;
  for(int r = 0; r < _rows; r++)
  {
    sum += _cell[r][col];
  }
  return sum;
}

int Matrix::cols() const
{
  return _cols;
}

int Matrix::rows() const
{
  return _rows;
}

void Matrix::__check(int row, int col) throw(EntropyException)
{
  if(row >= _rows || col >= _cols || row < 0 || col < 0)
  {
    std::stringstream oss;
    oss << "Invalid matrix index: ";
    if(row >= _rows)
    {
      oss << "row too large: " << row << " >= " << _rows;
    }
    if(row < 0)
    {
      oss << "row is negative: " << row;
    }
    if(col >= _cols)
    {
      oss << "col too large: " << col << " >= " << _cols;
    }
    if(col < 0)
    {
      oss << "col is negative: " << col;
    }
    throw EntropyException(oss.str());
  }
}

Matrix& Matrix::operator+=(const Matrix &m) throw(EntropyException)
{
  if(_rows != m.rows() || _cols != m.cols())
  {
    std::stringstream oss;
    oss << "Cannot += matrix, because ";
    if(_rows != m.rows())
    {
      oss << "number of rows do not match: " << _rows << " != " << m.rows();
    }
    if(_cols != m.cols())
    {
      oss << "number of cols do not match: " << _cols << " != " << m.cols();
    }
    throw EntropyException(oss.str());
  }

  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      __set(r, c, get(r,c) + m(r,c));
    }
  }
  return *this;
}

Matrix& Matrix::operator-=(const Matrix &m) throw(EntropyException)
{
  if(_rows != m.rows() || _cols != m.cols())
  {
    std::stringstream oss;
    oss << "Cannot -= matrix, because ";
    if(_rows != m.rows())
    {
      oss << "number of rows do not match: " << _rows << " != " << m.rows();
    }
    if(_cols != m.cols())
    {
      oss << "number of cols do not match: " << _cols << " != " << m.cols();
    }
    throw EntropyException(oss.str());
  }

  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      __set(r, c, get(r,c) - m(r,c));
    }
  }
  return *this;
}

Matrix& Matrix::operator*=(const double factor)
{
  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      __set(r, c, get(r, c) * factor);
    }
  }
  return *this;
}

const Matrix Matrix::operator+(const Matrix &m)
{
  Matrix r = *this;
  r += m;
  return r;
}

const Matrix Matrix::operator-(const Matrix &m)
{
  Matrix r = *this;
  r -= m;
  return r;
}

const Matrix Matrix::operator*(const double factor)
{
  Matrix r = *this;
  r *= factor;
  return r;
}


Matrix& Matrix::operator=(const Matrix &m)
{
  if(this == &m) return *this;

  if(_cell != NULL)
  {
    __deleteCells();
  }

  __init(m.rows(), m.cols(), 0.0);
  __copy(m);
  return *this;
}

void Matrix::__deleteCells()
{
  if (_cell != NULL)
  {
    for(register int i=_rows-1; i>=0; i--)
    {
      delete [] _cell[i];
    }
    delete [] _cell;
  }
}

void Matrix::__copy(const Matrix &m)
{
  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      __set(r, c, m(r,c));
    }
  }
}

void Matrix::setDiagonalMatrix(double value)
{
  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      if(r == c)
      {
        __set(r, c, value);
      }
      else
      {
        __set(r, c, 0.0);
      }
    }
  }
}

void Matrix::reset(int rows, int cols, double value)
{
  __deleteCells();
  __init(rows, cols, value);
}

Matrix& Matrix::operator=(const double d)
{
  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      __set(r, c, d);
    }
  }
  return *this;
}

Matrix& Matrix::operator-=(const double d)
{
  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      __set(r, c, get(r,c) - d);
    }
  }
  return *this;
}

Matrix& Matrix::operator/=(const double d)
{
  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      __set(r, c, get(r,c) / d);
    }
  }
  return *this;
}


void Matrix::rescaleRows(double value, bool verbose)
{
  for(int r = 0; r < _rows; r++)
  {
    double sum = 0;
    for(int c = 0; c < _cols; c++)
    {
      sum += get(r,c);
    }
    if(sum != value)
    {
      if(verbose)
      {
        std::cout << "row sum of row " << r << " is " << sum << std::endl;
      }
      sum /= value;
      for(int c = 0; c < _cols; c++)
      {
        __set(r, c, get(r,c) / sum);
      }
      if(verbose)
      {
        sum = 0;
        for(int c = 0; c < _cols; c++)
        {
          sum += get(r,c);
        }
        std::cout << "row sum of row " << r << " is now " << sum << std::endl;
      }
    }
  }
}

double Matrix::L2()
{
  double d = 0.0;
  double v = 0.0;
  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      v = get(r,c);
      d += v * v;
    }
  }
  return sqrt(d);
}

double Matrix::det() throw(EntropyException)
{
  if(cols() != rows())
  {
    throw EntropyException("Determinant is only implemented for quadratic matrices");
  }

  double result = 0; 

  if(_rows == 1)
  {
    result = get(0,0);
    return result; 
  } 

  if(_rows == 2)
  { 
    result = get(0, 0) * get(1, 1) - get(0, 1) * get(1, 0); 
    return result; 
  } 

  for(int i = 0; i < cols(); i++)
  { 
    Matrix tmp(_rows-1, _cols-1, 0.0);
    for(int j = 1; j < rows(); j++)
    { 
      for(int k = 0; k < cols(); k++)
      { 
        if(k < i)
        { 
          tmp(j - 1, k) = get(j,k); 
        }
        else if(k > i)
        { 
          tmp(j - 1, k - 1) = get(j,k); 
        } 
      } 
    } 
    result += get(0,i) * pow(-1, (double)i) * tmp.det();
  } 

  return result; 
}

void Matrix::transpose()
{
  Matrix m = (*this);
  reset(m.cols(), m.rows());

  for(int r = 0; r < m.rows(); r++)
  {
    for(int c = 0; c < m.cols(); c++)
    {
      __set(c,r, m(r,c));
    }
  }
}

void Matrix::adjunct()
{
  Matrix m(_rows, _cols);
  
  m = *this;

  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      Matrix ij = m;
      ij.cut(r,c);
      double v = pow(-1, (r+1)+(c+1)) * ij.det();
      __set(r, c, v);
    }
  }

  transpose();

}

void Matrix::cut(int r_index, int c_index)
{
  int r_new = _rows;
  int c_new = _cols;
  if(r_index >= 0)
  {
    r_new--;
  }
  if(c_index >= 0)
  {
    c_new--;
  }

  if(r_index < 0)
  {
    r_index = _rows + 1;
  }

  if(c_index < 0)
  {
    c_index = _cols + 1;
  }

  Matrix m(r_new, c_new);

  int r_offset = 0;
  int c_offset = 0;

  for(int i = 0; i < r_new; i++)
  {
    if(i < r_index)
    {
      r_offset = 0;
    }
    else
    {
      r_offset = 1;
    }
    for(int j = 0; j < c_new; j++)
    {
      if(j < c_index)
      {
        c_offset = 0;
      }
      else
      {
        c_offset = 1;
      }

      m(i, j) = get(i + r_offset, j + c_offset);
    }
  }

  *this = m;
}

Matrix& Matrix::operator*=(const Matrix &m)
{
  Matrix C(_rows, _cols);
  for(int r = 0; r < _rows; r++)
  {
    for(int c = 0; c < _cols; c++)
    {
      double v = 0.0;
      for(int k = 0; k < _rows; k++)
      {
        v += (get(r,k) * m(k,c));
      }
      C(r,c) = v;
    }
  }
  *this = C;
  return *this;
}

Matrix operator*(double factor, const Matrix& m) 
{ 
  Matrix R(m.rows(), m.cols());

  for(int r = 0; r < m.rows(); r++)
  {
    for(int c = 0; c < m.cols(); c++)
    {
      R(r, c) = m(r,c) * factor;
    }
  }
  return R;	
}

// from http://users.erols.com/mdinolfo/matrix.htm
void Matrix::invert()
{
  if(_rows != _cols && _cols <=1)
  {
    throw EntropyException("Matrix is not rectangular or too small");
  }

  for (int i = 1; i < _rows; i++ )
  {
    __set(0, i, __get(0,i) / __get(0, 0));
  }
  for (int i = 1; i < _rows; i++ )
  {
    for (int j = i; j < _rows; j++ ) // do a column of L
    {
      double sum = 0.0;
      for (int k = 0; k < i; k++ )
      {
        sum += __get(j, k) * __get(k,i);
      }
      __add(j,i,-sum);
    }
    if ( i == _rows - 1 ) continue;
    for (int j = i + 1; j < _rows; j++ )  // do a row of U
    {
      double sum = 0.0;
      for (int k = 0; k < i; k++ )
      {
        sum += __get(i,k) * __get(k,j);
      }
      __set(i, j, (__get(i,j) - sum ) / __get(i,i));
    }
  }
  for (int i = 0; i < _rows; i++ ) // invert L
  {
    for (int j = i; j < _rows; j++ )
    {
      double x = 1.0;
      if ( i != j )
      {
        x = 0.0;
        for (int k = i; k < j; k++ )
        {
          x -= __get(j,k) * __get(k,i);
        }
      }
      __set(j,i,  x / __get(j,j));
    }
  }
  for (int i = 0; i < _rows; i++ )  // invert U
  {
    for (int j = i; j < _rows; j++ )
    {
      if ( i == j ) continue;
      double sum = 0.0;
      for (int k = i; k < j; k++ )
      {
        sum += __get(k,j) * ((i == k)?1.0:__get(i,k));
      }
      __set(i,j, -sum);
    }
  }
  for (int i = 0; i < _rows; i++ )  // final inversion
  {
    for (int j = 0; j < _rows; j++ )
    {
      double sum = 0.0;
      for (int k = ( ( i > j ) ? i : j ); k < _rows; k++ )
      {
        sum += ( ( j == k ) ? 1.0 : __get(j,k)) * __get ( k, i );
      }
      __set(j, i, sum);
    }
  }
}

double Matrix::__get(const int row, const int col)
{
  return _cell[row][col];
}

void Matrix::__add(const int row, const int col, const double value)
{
  _cell[row][col] += value;
}

const Matrix Matrix::operator*(const Matrix &m)
{
  Matrix r = *this;
  r *= m;
  return r;
}

