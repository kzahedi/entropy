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

      static DContainer* read(string filename, int n, ...);
      static DContainer* read(string filename);
      static DContainer* read(string filename, vector<int>);
      static void write(string filename, DContainer* c);
      static void write(string filename, ULContainer* c);
      static void write(string filename, IContainer* c);

    private:
  };
}

#endif // __CSV_H__
