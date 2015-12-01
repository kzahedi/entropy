#include <entropy++/CsvToContainer.h>

#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>

using namespace boost;

CsvToContainer::CsvToContainer()
{
}

Container* CsvToContainer::read(string filename, int n, ...)
{

  vector<int> indices;
  va_list ap;
  va_start(ap, n);
  for(int i = 0; i < n; i++)
  {
     indices.push_back(va_arg(ap, int));
  }
  va_end(ap);

  ifstream ifs(filename.c_str());
  string   line;
  int      nrOfLines = 0;
  while (std::getline(ifs, line))
  {
    nrOfLines++;
  }
  ifs.close();
  ifs.open(filename.c_str());

  Container *c = new Container(nrOfLines, indices.size());

  while (std::getline(ifs, line))
  {
    boost::char_separator<char> sep(",");
    tokenizer<boost::char_separator<char> > tk(line,sep);
    vector<string> vec;
    for(tokenizer<boost::char_separator<char> >::iterator i(tk.begin()); i!=tk.end();++i) 
    {
      vec.push_back(*i);
    }

    for(vector<int>::iterator i = indices.begin(); i != indices.end(); i++)
    {
      (*c) << atof(vec[*i].c_str());
    }
  }

  return c;
}
