#include <entropy++/Csv.h>

#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <cstdarg>

using namespace boost;

Csv::Csv()
{
}

DContainer* Csv::read(string filename)
{
  ifstream ifs(filename.c_str());
  string   line;
  int      nrOfLines   = 1;
  int      nrOfColumns = -1;
  if(std::getline(ifs, line))
  {
    boost::char_separator<char> sep(",");
    tokenizer<boost::char_separator<char> > tk(line,sep);
    vector<string> vec;
    for(tokenizer<boost::char_separator<char> >::iterator i(tk.begin()); i!=tk.end();++i) 
    {
      vec.push_back(*i);
    }
    nrOfColumns = (int)vec.size();
  }

  while (std::getline(ifs, line))
  {
    nrOfLines++;
  }
  ifs.close();
  ifs.open(filename.c_str());

  DContainer *c = new DContainer(nrOfLines, nrOfColumns);

  while (std::getline(ifs, line))
  {
    boost::char_separator<char> sep(",");
    tokenizer<boost::char_separator<char> > tk(line,sep);
    for(tokenizer<boost::char_separator<char> >::iterator i(tk.begin()); i!=tk.end();++i) 
    {
      (*c) << atof((*i).c_str());
    }
  }

  return c;
}

DContainer* Csv::read(string filename, int n, ...)
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

  DContainer *c = new DContainer(nrOfLines, indices.size());

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

DContainer* Csv::read(string filename, vector<int> indices)
{
  ifstream ifs(filename.c_str());
  string   line;
  int      nrOfLines = 0;
  while (std::getline(ifs, line))
  {
    nrOfLines++;
  }
  ifs.close();
  ifs.open(filename.c_str());

  DContainer *c = new DContainer(nrOfLines, indices.size());

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
void Csv::write(string filename, DContainer* container)
{
  ofstream out(filename.c_str());
  for(int r = 0; r < container->rows(); r++)
  {
    out << (*container)(r,0);
    for(int c = 1; c < container->columns(); c++)
    { 
      out << "," << (*container)(r,c);
    }
    out << endl;
  }
  out.close();
}

void Csv::write(string filename, ULContainer* container)
{
  ofstream out(filename.c_str());
  for(int r = 0; r < container->rows(); r++)
  {
    out << (*container)(r,0);
    for(int c = 1; c < container->columns(); c++)
    { 
      out << "," << (*container)(r,c);
    }
    out << endl;
  }
  out.close();
}
