#ifndef __CSV_TO_CONTAINER_H__
#define __CSV_TO_CONTAINER_H__

#include <entropy++/Container.h>

#include <vector>
#include <string>

using namespace std;

class CsvToContainer
{
  public:
    CsvToContainer();

    Container* read(string filename, int n, ...);

  private:
};

#endif // __CSV_TO_CONTAINER_H__
