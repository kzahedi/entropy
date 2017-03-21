#ifndef __SPARSE_MATRIX_H__
#define __SPARSE_MATRIX_H__

#include <entropy++/EntropyException.h>

#include <vector>
#include <string>
#include <iostream>
#include <utility>

using namespace std;

namespace entropy
{
  class MatrixIndex
  {
    public:
      MatrixIndex();
      // ~MatrixIndex();

      MatrixIndex(int f)
      {
        first  =  f;
        second = -1;
        third  = -1;
      }

      MatrixIndex(int f, int s)
      {
        first  =  f;
        second =  s;
        third  = -1;
      }

      MatrixIndex(int f, int s, int t)
      {
        first  = f;
        second = s;
        third  = t;
      }

      MatrixIndex(const MatrixIndex& m)
      {
        first  = m.first;
        second = m.second;
        third  = m.third;
      }

      MatrixIndex operator=(const MatrixIndex &m)
      {
        this->first  = m.first;
        this->second = m.second;
        this->third  = m.third;
        return *this;
      }

      int first;
      int second;
      int third;

    private:
  };


  class SparseMatrix
  {
    public:

      /** Most simple constructor. Every cell is initialised with 0.*/
      SparseMatrix();
      SparseMatrix(double);
      SparseMatrix(SparseMatrix &m);
      SparseMatrix(const SparseMatrix &m);

      /** Destructor. */
      ~SparseMatrix();

      /** Allows to access and modification of the values by indexing */
      double&   operator()(int first) throw(EntropyException);
      double    operator()(int first) const throw(EntropyException);

      double&   operator()(int first, int second) throw(EntropyException);
      double    operator()(int first, int second) const throw(EntropyException);

      double&   operator()(int first, int second, int third) throw(EntropyException);
      double    operator()(int first, int second, int third) const throw(EntropyException);

      SparseMatrix&   operator+=(const SparseMatrix &m)  throw(EntropyException);
      SparseMatrix&   operator-=(const SparseMatrix &m)  throw(EntropyException);

      // TODO: SparseMatrix A(10,10); SparseMatrix A = B; does not work
      const SparseMatrix operator* (const double factor);
      SparseMatrix&      operator*=(const double factor);

      SparseMatrix&      operator=(const SparseMatrix &m);

      SparseMatrix&      operator-=(const double d);
      SparseMatrix&      operator/=(const double d);

      const SparseMatrix operator-(const SparseMatrix &m);
      const SparseMatrix operator+(const SparseMatrix &m);

      int size() const;
      double get(int);
      double value(int,int,int);
      double value(int,int);
      double value(int);
      double defaultValue();
      bool available(int first, int second, int third);
      bool available(int first, int second);
      bool available(int first);

      MatrixIndex getmi(int) const;

      void reset();

      // double L2();

      // int    cols() const;
      // int    rows() const;

      // void   setDiagonalSparseMatrix(double value);

      // double rowSum(const int row);
      // double colSum(const int col);

      // double l2();

      // void   rescaleRows(double value, bool verbose);

      // double det() throw(EntropyException);
      // void   invert();
      // void   transpose();
      // void   adjunct();
      // void   cut(int r_index = -1, int c_index = -1);

      // friend SparseMatrix operator*(double, const SparseMatrix&);

      // friend std::ostream& operator<<(std::ostream& str, const SparseMatrix& m)
      // {
      // for(int r = 0; r < m._rows; r++)
      // {
      // for(int c = 0; c < m._cols - 1; c++)
      // {
      // str << m._cell[r][c] << " ";
      // }
      // str << m._cell[r][m._cols -1] << std::endl;
      // }
      // return str;
      // };

    protected:
      void   __set(const int f, const int s, const int t, const double value);
      double __get(const int f, const int s, const int t) const;

      void   __copy(const SparseMatrix &m);

      vector<MatrixIndex> _indices;
      vector<double>      _values;

      int    _size;
      double _default;

  };
}
#endif // __SPARSE_MATRIX_H__
