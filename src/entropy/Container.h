#ifndef __CONTAINER_H__
#define __CONTAINER_H__

# define CONTAINER_DISCRETISE_UNIFORM 1001

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

    Container* discretise();

  private:
    Container* __uniformDiscretisation();
    int        __discretiseAndCombineValues(double *values);

    double** _data;
    int      _mode;
    int      _rows;
    int      _columns;
    int      _fillIndex;
    int*     _bins;
    double** _domains;

    bool     _domainsGiven;
    bool     _binsGiven;
};

#endif // __CONTAINER_H__
