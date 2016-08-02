#ifndef __CSV_H__
#define __CSV_H__

#include <entropy++/Container.h>

#include <vector>
#include <string>

using namespace std;

class Csv
{
  public:
    Csv();

    Container* read(string filename, int n, ...);
    Container* read(string filename);
    Container* read(string filename, vector<int>);
    void write(string filename, Container* c);

  private:
};

#endif // __CSV_H__
