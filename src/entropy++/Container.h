#ifndef __CONTAINER_H__
#define __CONTAINER_H__

#include <entropy++/defs.h>

#include <assert.h>
#include <cstdarg>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <ostream>
#include <vector>
#include <algorithm>

# define MIN(a,b) (((a)<(b))?a:b)
# define MAX(a,b) (((a)>(b))?a:b)

using namespace std;

#define FILL_MODE_BY_ROW    1001
#define FILL_MODE_BY_COLUMN 1002

using namespace std;

namespace entropy
{
  template <class T>
    class Container 
    {
      public:

        Container<T>(int rows, int columns)
        {
          _data         = new T*[rows];
          _domains      = new double*[columns];
          _bins         = new int[columns];
          _mode         = UNIFORM;
          _fillMode     = FILL_MODE_BY_ROW;

          _domainsGiven = false;
          _binsGiven    = false;
          _discretised  = false;

          for(int i = 0; i < columns; i++)
          {
            _domains[i] = new double[2];
          }

          for(int i = 0; i < rows; i++)
          {
            _data[i] = new T[columns];
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

        ~Container()
        {
          if(_data != NULL)
          {
            for(int r = 0; r < _rows; r++) delete _data[r];
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


        // Container(const Container);
        Container<T>& operator=(const Container<T>& c)
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

        const Container& operator<<(const T& value) const
        {
          assert(_fillIndex < _rows * _columns);

          int c = -1;
          int r = -1;

          switch(_fillMode)
          {
            case FILL_MODE_BY_ROW:
              c = _fillIndex % _columns;
              r = _fillIndex / _columns;
              break;
            case FILL_MODE_BY_COLUMN:
              c = _fillIndex / _rows;
              r = _fillIndex % _rows;
              break;
          }
          _data[r][c] = value;
          Container *co = (Container*)this;
          co->_fillIndex++;

          return *this;
        }

        // merge
        Container<T>& operator+=(const Container<T>& c)
        {
          int newrows    = MIN(this->_rows, c._rows);
          int newcolumns = this->_columns + c._columns;

          T** tmp      = new T*[newrows];
          int* newBins = new int[newcolumns];

          for(int i = 0; i < _columns; i++)
          {
            newBins[i] = _bins[i];
          }
          for(int i = 0; i < c._columns; i++)
          {
            newBins[_columns + i] = c._bins[i];
          }

          for(int i = 0; i < newrows; i++)
          {
            tmp[i] = new T[newcolumns];
            int index = 0;
            for(int j = 0; j < _columns; j++)
            {
              tmp[i][index] = _data[i][j];
              index++;
            }
            for(int j = 0; j < c._columns; j++)
            {
              tmp[i][index] = c(i,j);
              index++;
            }
          }

          if(_data != NULL)
          {
            for(int r = 0; r < _rows; r++) delete _data[r];
            delete _data;
          }

          _data    = tmp;
          _columns = newcolumns;

          delete _bins;
          _bins = newBins;

          return *this;

        }

        T  operator()(const int row, const int column) const
        {
          // assert(row    < _rows);
          // assert(column < _columns);
          return _data[row][column];
        }

        T& operator()(const int row, const int column)
        {
          // assert(row    < _rows);
          // assert(column < _columns);
          return _data[row][column];
        }

        T get(int row, int column)
        {
          // assert(row    < _rows);
          // assert(column < _columns);
          return _data[row][column];
        }

        void set(int row, int column, T value)
        {
          // assert(row    < _rows);
          // assert(column < _columns);
          _data[row][column] = value;
        }

        void normaliseColumn(int c, double min, double max)
        {
          for(int r = 0; r < _rows; r++)
          {
            _data[r][c] = (_data[r][c] - min) / (max - min);
          }
        }

        void setFillMode(int f)
        {
          _fillMode = f;
        }

        T rowSum(int r)
        {
          T sum = (T)0;
          for(int c = 0; c < _columns; c++)
          {
            sum += _data[r][c];
          }
          return sum;
        }

        T colSum(int c)
        {
          T sum = (T)0;
          for(int r = 0; r < _rows; r++)
          {
            sum += _data[r][c];
          }
          return sum;
        }

        double max()
        {
          double m = max(0);
          for(int i = 1; i < _columns; i++)
          {
            double n = max(i);
            if(n > m) m = n;
          }
          return m;
        }

        double max(int column)
        {
          double m = _data[0][column];
          for(int i = 1; i < _rows; i++)
          {
            if(_data[i][column] > m) m = _data[i][column];
          }
          return m;
        }

        double min()
        {
          double m = min(0);
          for(int i = 1; i < _columns; i++)
          {
            double n = min(i);
            if(n < m) m = n;
          }
          return m;
        }

        double min(int column)
        {
          double m = _data[0][column];
          for(int i = 1; i < _rows; i++)
          {
            if(_data[i][column] < m) m = _data[i][column];
          }
          return m;
        }

        int rows()
        {
          return _rows;
        }

        int columns()
        {
          return _columns;
        }

        Container<T>* columns(int n, ...)
        {
          vector<int> indices;
          va_list ap;
          va_start(ap, n);
          for(int i = 0; i < n; i++)
          {
            indices.push_back(va_arg(ap, int));
          }
          va_end(ap);

          Container<T> *extracted = new Container<T>(this->rows(), n);

          for(int r = 0; r < _rows; r++)
          {
            for(int i = 0; i < (int)indices.size(); i++)
            {
              (*extracted) << get(r, indices[i]);
            }
          }

          return extracted;
        }

        Container<T>* copy()
        {
          Container<T> *copy = new Container<T>(_rows, _columns);
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

        bool isDiscretised()
        {
          return _discretised;
        }

        void isDiscretised(bool b)
        {
          _discretised = b;
        }

        Container<T>* drop(int n)
        {
          if(n > 0) return __dropFirst(n);
          if(n < 0) return __dropLast(n);
          return copy();
        }

        void setBinSizes(int bins)
        {
          for(int c = 0; c < _columns; c++)
          {
            _bins[c] = bins;
          }
          _binsGiven = true;
        }

        void setBinSizes(int *bins)
        {
          for(int c = 0; c < _columns; c++)
          {
            _bins[c] = bins[c];
          }
          _binsGiven = true;
        }

        void setDomains(double **domains)
        {
          for(int c = 0; c < _columns; c++)
          {
            _domains[c][0] = domains[c][0];
            _domains[c][1] = domains[c][1];
          }
          _domainsGiven = true;
        }

        void setDomains(double min, double max)
        {
          for(int c = 0; c < _columns; c++)
          {
            _domains[c][0] = min;
            _domains[c][1] = max;
          }
          _domainsGiven = true;
        }

        void setDiscretisationMode(int mode)
        {
          _mode = mode;
        }

        T columnSum(int c)
        {
          T s = 0.0;
          for(int r = 0; r < _rows; r++) s += _data[r][c];
          return s;
        }

        Container<unsigned long>* discretise()
        {
          Container<unsigned long> *dis = NULL;
          Container<unsigned long> *com = NULL;
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

        Container<unsigned long>* discretiseByColumn()
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

        Container<unsigned long>* combineDiscretisedColumns()
        {
          assert(_discretised == true);
          Container<unsigned long> *copy = new Container<unsigned long>(_rows, 1);
          // __copyProperties(copy);
          copy->isDiscretised(_discretised);

          int maxBins[_columns];
          for(int c = 0; c < _columns; c++)
          {
            maxBins[c] = max(c);
          }
          for(int r = 0; r < _rows; r++)
          {
            unsigned long v = 0;
            unsigned long f = 1;
            for(int c = 0; c < _columns; c++)
            {
              if(c > 0) f = f * MAX(maxBins[c-1],1);
              v += f * get(r,c);
            }
            *copy << v;
          }
          copy->relabel();
          return copy;
        }

        friend std::ostream& operator<<(std::ostream& str, const Container<T>& container)
        {
          str << "Container content:" << endl;
          for(int r = 0; r < container._rows; r++)
          {
            str << "  " << container(r,0); 
            for(int c = 1; c < container._columns; c++)
            {
              str << " " << container(r, c);
            }
            str << endl;
          }
          return str;
        };

        void relabel()
        {
          for(int c = 0; c < _columns; c++)
          {
            vector<int> values;
            for(int r = 0; r < _rows; r++)
            {
              int  value = (int)_data[r][c];
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
              int value = (int)_data[r][c];
              for(int i = 0; i < (int)values.size(); i++)
              {
                if(value == values[i])
                {
                  _data[r][c] = (double)i;
                }
              }
            }
          }
        }

        // returns a single column container with the unique values taken over all
        // columns
        Container<T>* globalUnique()
        {
          // Container<T> *new = new Container<T>(this->rows(), 1);
          vector<T> values;
          for(int c = 0; c < this->columns(); c++)
          {
            for(int r = 0; r < this->rows(); r++)
            {
              T value = this->get(r,c);
              if(std::find(values.begin(), values.end(), value) == values.end())
                //if(values.find(value) == values.end())
              {
                values.push_back(value);
              }
            }
          }
          Container<T>* new_container = new Container<T>(values.size(), 1);
          for(int i = 0; i < values.size(); i++)
          {
            *new_container << values[i];
          }
          return new_container;
        }

        Container<T>* unique()
        {
          // Container<T> *new = new Container<T>(this->rows(), 1);
          vector<int> uniqueIndices;
          for(int row1 = 0; row1 < this->rows(); row1++)
          {
            bool found = false;
            for(int row2 = row1 + 1; row2 < this->rows(); row2++)
            {
              bool rowEqual = true;
              for(int c = 0; c < this->columns(); c++)
              {
                rowEqual &= (this->get(row1, c) == this->get(row2, c));
                if(rowEqual == false) break;
              }
              if(rowEqual == true)
              {
                found = true;
                break;
              }
            }
            if(found == false) uniqueIndices.push_back(row1);
          }

          Container<T>* new_container = new Container<T>(uniqueIndices.size(), this->columns());
          for(int i = 0; i < uniqueIndices.size(); i++)
          {
            for(int c = 0; c < this->columns(); c++)
            {
              *new_container << this->get(uniqueIndices[i], c);
            }
          }
          return new_container;
        }

        int find(T* values)
        {
          bool found = false;
          for(int r = 0; r < this->rows(); r++)
          {
            found = true;
            for(int c = 0; c < this->columns(); c++)
            {
              found &= (this->get(r,c) == values[c]);
              if(found == false) break;
            }
            if(found == true) return r;
          }
          return -1;
        }

        int find(Container<T>* other, int row)
        {
          bool found = false;
          for(int r = 0; r < this->rows(); r++)
          {
            found = true;
            for(int c = 0; c < this->columns(); c++)
            {
              found &= (this->get(r,c) == other->get(row,c));
              if(found == false) break;
            }
            if(found == true) return r;
          }
          return -1;
        }

        vector<int> findlist(T* values)
        {
          vector<int> indices;
          bool found = false;
          for(int r = 0; r < this->rows(); r++)
          {
            found = true;
            for(int c = 0; c < this->columns(); c++)
            {
              found &= (this->get(r,c) == values[c]);
              if(found == false) break;
            }
            if(found == true) indices.push_back(r);
          }
          return indices;
        }

        vector<int> findlist(Container<T>* other, int row)
        {
          vector<int> indices;
          bool found = false;
          for(int r = 0; r < this->rows(); r++)
          {
            found = true;
            for(int c = 0; c < this->columns(); c++)
            {
              found &= (this->get(r,c) == other->get(row,c));
              if(found == false) break;
            }
            if(found == true) indices.push_back(r);
          }
          return indices;
        }

      private:
        Container<unsigned long>* __uniformDiscretisationByColumn()
        {
          Container<unsigned long> *copy = new Container<unsigned long>(_rows, _columns);
          // __copyProperties(copy);
          copy->isDiscretised(_discretised);
          for(int c = 0; c < _columns; c++)
          {
            vector<int> values;
            for(int r = 0; r < _rows; r++)
            {
              ASSERT(_domains[c][0] <= get(r,c) && get(r,c) <= _domains[c][1], "get(" << r << "," << c << ") = " << get(r,c) << endl << "Domain " << _domains[c][0] << ", " << _domains[c][1] << endl);
              unsigned long mapped  = (unsigned long)(((get(r,c)       - _domains[c][0])
                                                       / (_domains[c][1] - _domains[c][0]))
                                                      * _bins[c]);
              unsigned long cropped = (unsigned long)MIN(_bins[c]-1, mapped);
              copy->set(r, c, cropped);
            }
          }
          copy->relabel();
          copy->isDiscretised(true);
          return copy;
        }


        void __copyProperties(Container* dst)
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

        Container<T>* __dropFirst(int n)
        {
          Container<T> *copy = new Container<T>(_rows-abs(n), _columns);
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

        Container<T>* __dropLast(int n)
        {
          Container<T> *copy = new Container<T>(_rows-abs(n), _columns);
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

        T**         _data;
        int         _mode;
        int         _rows;
        int         _columns;
        int         _fillIndex;
        int*        _bins;
        int         _fillMode;
        double**    _domains;


        bool        _domainsGiven;
        bool        _binsGiven;
        bool        _discretised;
    };

  typedef Container<float> FContainer;
  typedef Container<double> DContainer;
  typedef Container<int> IContainer;
  typedef Container<unsigned long> ULContainer;
}

#endif // __CONTAINER_H__
