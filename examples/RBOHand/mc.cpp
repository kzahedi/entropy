#include <iostream>
#include <fstream>

#include <stdlib.h>

#include <entropy++/Csv.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/sparse/state/MC_W.h>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <string>

#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>

#define SET_DOMAIN_AND_BINS(a,d,b) \
  a->setDomains(d);   \
  a->setBinSizes(b);

using namespace std;
using namespace boost;

double x_min;
double x_max;

double y_min;
double y_max;

double z_min;
double z_max;

vector<double> amin;
vector<double> amax;

vector<string> files;

DEFINE_string(d,    "",   "directory");
DEFINE_string(wi,   "",   "W column indices");
DEFINE_string(ai,   "",   "A column indices");
DEFINE_int64(wbins, 1000, "world state bins");
DEFINE_int64(abins, 100,  "action bins");

vector<int> int_tokenizer(string s)
{
  boost::char_separator<char> sep(",");
  tokenizer<boost::char_separator<char> > tk(s,sep);
  vector<int> vec;
  for(tokenizer<boost::char_separator<char> >::iterator i(tk.begin()); i!=tk.end();++i) 
  {
    vec.push_back(atoi((*i).c_str()));
  }
  return vec;
}

int main(int argc, char** argv)
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;

  VLOG(10) << "Directory: " << FLAGS_d;

  string w_domain_file = string(FLAGS_d) + "/W.domains.csv";
  string a_domain_file = string(FLAGS_d) + "/A.domains.csv";
  string w_states      = string(FLAGS_d) + "/hand.sofastates.csv";
  string a_states      = string(FLAGS_d) + "/control.states.csv";

  cout << "Files:" << endl;
  cout << "  W domain: " << w_domain_file << endl;
  cout << "  A domain: " << a_domain_file << endl;
  cout << "  W states: " << w_states      << endl;
  cout << "  A states: " << a_states      << endl;

  Csv* csv = new Csv();

  Container* w_domains = NULL;
  if(FLAGS_wi == "")
  {
    w_domains = csv->read(w_domain_file);
  }
  else
  {
    vector<int> w_indices = int_tokenizer(FLAGS_wi);
    w_domains = csv->read(w_domain_file, w_indices);
  }

  Container* a_domains = NULL;
  if(FLAGS_ai == "")
  {
    a_domains = csv->read(a_domain_file);
  }
  else
  {
    vector<int> a_indices = int_tokenizer(FLAGS_ai);
    a_domains = csv->read(a_domain_file, a_indices);
  }

  cout << "W domain file:" << endl;
  cout << *w_domains << endl;

  cout << "A domain file:" << endl;
  cout << *a_domains << endl;

  // double** w_domains = new double[e->wsize()];
  // for(int i = 0; i < e->wsize(); i++)
  // {
    // w_domains[i] = new double[2];
    // w_domains[i][0] = e->wmin(0,i);
    // w_domains[i][1] = e->wmax(0,i);
  // }




  // Container* W1 =  W->drop(-1);
  // Container* A1 = dA->drop(-1);
  // Container* W2 =  W->drop(W->rows() - A1->rows());

  // double     mc  = entropy::sparse::MC_W(W2, W1, A1);
  // Container* mcd = entropy::sparse::state::MC_W(W2, W1, A1);

  // sst.str("");
  // sst << FLAGS_out << "-state-dependent.csv";
  // csv->write(sst.str().c_str(), mcd);
  // delete csv;

  // cout << "MC: " << mc << endl;

  // sst.str("");
  // sst << FLAGS_out << "-averaged.txt";
  // VLOG(10) << "writing \"" << sst.str() << "\"";
  // ofstream output;
  // output.open(sst.str());
  // output << mc << endl;
  // output.close();
}
