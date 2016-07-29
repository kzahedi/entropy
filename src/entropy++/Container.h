#ifndef __CONTAINER_H__
#define __CONTAINER_H__

#include <entropy++/defs.h>

#include <ostream>
#include <vector>

#define FILL_MODE_BY_ROW    1001
#define FILL_MODE_BY_COLUMN 1002

using namespace std;

class Container 
{
  public:

    Container(int rows, int columns);
    ~Container();

    // Container(const Container);
    Container& operator=(const Container&);

    const Container& operator<<(const double&) const;
    // merge
    Container& operator+=(const Container&);

    double  operator()(const int row, const int column) const;
    double& operator()(const int row, const int column);

    double     get(int, int);
    void       set(int, int, double);
    void       normaliseColumn(int, double, double);

    void       setFillMode(int);

    double     max();
    double     max(int);

    double     min();
    double     min(int);

    int        rows();
    int        columns();
    Container* columns(int n, ...);
    Container* copy();

    bool       isDiscretised();
    void       isDiscretised(bool);
    Container* drop(int n);
    // void take(int n);

    void setBinSizes(int*);
    void setDomains(double**);
    void setDiscretisationMode(int);

    double columnSum(int);

    Container* discretise();
    Container* discretiseByColumn();
    Container* combineDiscretisedColumns();

    friend std::ostream& operator<<(std::ostream& str, const Container& container)
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

  private:
    Container* __uniformDiscretisation();
    Container* __uniformDiscretisationByColumn();
    int        __discretiseAndCombineValues(double *values);
    void       __strip();
    void       __copyProperties(Container* dst);
    Container* __dropFirst(int n);
    Container* __dropLast(int n);

    double**    _data;
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

#endif // __CONTAINER_H__
