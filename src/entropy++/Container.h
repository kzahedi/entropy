#ifndef __CONTAINER_H__
#define __CONTAINER_H__

# define CONTAINER_DISCRETISE_UNIFORM 1001

#include <ostream>

using namespace std;

class Container 
{
  public:

    Container(int rows, int columns);
    ~Container();

    const Container& operator<<(const double&) const;
    double operator()(const int row, const int column) const;
    double get(int, int);

    int rows();
    int columns();

    void setBinSizes(int*);
    void setDomains(double**);
    void setDiscretisationMode(int);

    bool isDiscretised();
    Container* discretise();

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
    int        __discretiseAndCombineValues(double *values);
    void       __strip();

    double** _data;
    int      _mode;
    int      _rows;
    int      _columns;
    int      _fillIndex;
    int*     _bins;
    double** _domains;

    bool     _domainsGiven;
    bool     _binsGiven;
    bool     _discretised;
};

#endif // __CONTAINER_H__
