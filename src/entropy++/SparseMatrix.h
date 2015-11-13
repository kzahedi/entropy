#ifndef __SPARSE_MATRIX_H__
#define __SPARSE_MATRIX_H__

#include <entropy++/MatrixException.h>

#include <vector>
#include <string>
#include <iostream>
#include <utility>

using namespace std;

typedef pair<int, int> matrixindex;

class SparseMatrix 
{
  public:
    /** Most simple constructor. Every cell is initialised with 0.*/
    SparseMatrix();
    SparseMatrix(SparseMatrix &m);
    SparseMatrix(const SparseMatrix &m);

    /** Destructor. */
    ~SparseMatrix();

    /** Allows to access and modification of the values by indexing */
    double&   operator()(int row, int col) throw(MatrixException);
    double    operator()(int row, int col) const throw(MatrixException);
    SparseMatrix&   operator+=(const SparseMatrix &m)  throw(MatrixException);
    SparseMatrix&   operator-=(const SparseMatrix &m)  throw(MatrixException);

    // TODO: SparseMatrix A(10,10); SparseMatrix A = B; does not work
    const SparseMatrix operator* (const double factor);
    SparseMatrix&      operator*=(const double factor);

    SparseMatrix&      operator= (const SparseMatrix &m);
    SparseMatrix&      operator*=(const SparseMatrix &m);

    SparseMatrix&      operator-=(const double d);
    SparseMatrix&      operator/=(const double d);

    const SparseMatrix operator-(const SparseMatrix &m);
    const SparseMatrix operator+(const SparseMatrix &m);

    double get(const int row, const int col);

    int size() const;

    matrixindex getmi(int) const;

    // double L2();

    // int    cols() const;
    // int    rows() const;

    // void   setDiagonalSparseMatrix(double value);

    // double rowSum(const int row);
    // double colSum(const int col);

    // double l2();

    void reset();
    // void   rescaleRows(double value, bool verbose);

    // double det() throw(MatrixException);
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
    void     __set(const  int row, const int col, const  double value);
    double   __get(const  int row, const int col) const;
    void     __copy(const SparseMatrix &m);
    int      _rows;
    int      _cols;
    double   _tmp;

    vector<matrixindex> _indices;
    vector<double>      _values;

};
#endif // __SPARSE_MATRIX_H__
