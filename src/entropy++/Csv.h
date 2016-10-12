#ifndef __CSV_H__
#define __CSV_H__

#include <entropy++/Container.h>

#include <vector>
#include <string>

using namespace std;

namespace entropy
{
  class Csv
  {
    public:
      Csv();

      DContainer* read(string filename, int n, ...);
      DContainer* read(string filename);
      DContainer* read(string filename, vector<int>);
      void write(string filename, DContainer* c);
      void write(string filename, ULContainer* c);

    private:
  };
}

#endif // __CSV_H__
