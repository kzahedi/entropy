#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <math.h>

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

DEFINE_string(d,    "",    "directory");
DEFINE_string(wi,   "",    "W column indices");
DEFINE_string(ai,   "",    "A column indices");
DEFINE_int64(wbins, 1000,  "world state bins");
DEFINE_int64(abins, 100,   "action bins");
DEFINE_bool(csv,    false, "write all data that is used for the analysis into csv files");

void check_domains(string label, DContainer *domain)
{
  bool found = false;
  for(int i = 0; i < domain->columns(); i++)
  {
    if(fabs((*domain)(0,i) - (*domain)(1,i)) < 0.00001)
    {
      found = true;
      cerr << "ERROR: In " << label << " domains: found " << 
        "[" << (*domain)(0,i) << ", " << (*domain)(1,i) << "] at index " << i << endl;
    }
  }
  if(found) 
  {
    cout << *domain << endl;
    exit(-1);
  }
}

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
  FLAGS_log_dir = ".";

  // if(FLAGS_lf == false)
  // {
    // FLAGS_logtostderr = 1;
  // }

  VLOG(10) << "Directory: " << FLAGS_d;

  string w_domain_file = string(FLAGS_d) + "/W.domains.csv";
  string a_domain_file = string(FLAGS_d) + "/A.domains.csv";
  string w_states      = string(FLAGS_d) + "/hand.sofastates.csv";
  string a_states      = string(FLAGS_d) + "/control.states.csv";

  VLOG(1) << "Files:";
  VLOG(1) << "  W domain: " << w_domain_file;
  VLOG(1) << "  A domain: " << a_domain_file;
  VLOG(1) << "  W states: " << w_states;
  VLOG(1) << "  A states: " << a_states;

  Csv* csv = new Csv();

  // W data
  DContainer* w_domains = NULL;
  DContainer* W         = NULL;
  double w_min[3];
  double w_max[3];
  vector<int> w_indices;
  if(FLAGS_wi == "")
  {
    w_domains = csv->read(w_domain_file);
    W         = csv->read(w_states);
    VLOG(100) << "W domain:";
    VLOG(100) << *w_domains;
    for(int i = 0; i < 3; i++)
    {
      w_min[i]  = (*w_domains)(0,i);
      w_max[i]  = (*w_domains)(1,i);
    }
    for(int i = 0; i < w_domains->columns(); i++)
    {
      if((*w_domains)(0,i) < w_min[i%3]) w_min[i%3] = (*w_domains)(0,i);
      if((*w_domains)(1,i) > w_max[i%3]) w_max[i%3] = (*w_domains)(1,i);
    }
  }
  else
  {
    w_indices = int_tokenizer(FLAGS_wi);
    w_domains = csv->read(w_domain_file, w_indices);
    W         = csv->read(w_states, w_indices);
    VLOG(100) << "W domain:";
    VLOG(100) << *w_domains;
    for(int i = 0; i < 3; i++)
    {
      w_min[i] = -1000000.0;
      w_max[i] = -1000000.0;
    }
    for(int i = 0; i < (int)w_indices.size(); i++)
    {
      // map the chosen indices to the x,y,z coordinates
      int index = w_indices[i] % 3;
      if(w_min[index] > (*w_domains)(0, i) ||
         w_min[index] == -1000000.0)
        w_min[index] = (*w_domains)(0, i);
      if(w_max[index] < (*w_domains)(1, i) ||
         w_max[index] == -1000000.0)
        w_max[index] = (*w_domains)(1, i);
    }
  }


  check_domains("W", w_domains);
  VLOG(1) << "W min/max:";
  for(int i = 0; i < 3; i++) VLOG(1) << w_min[i] << ", " << w_max[i];

  VLOG(1) << "W domain file:";
  VLOG(1) << *w_domains;
  VLOG(1) << "W states:";
  VLOG(1) << *W;

  double **wdomain = new double*[W->columns()];
  int *wbins = new int[W->columns()];
  for(int i = 0; i < W->columns(); i++) wbins[i] = FLAGS_wbins;

  if(FLAGS_wi == "")
  {
    for(int i = 0; i < W->columns(); i++)
    {
      wdomain[i]    = new double[2];
      wdomain[i][0] = w_min[i%3];
      wdomain[i][1] = w_max[i%3];
    }
  }
  else
  {
    for(int i = 0; i < W->columns(); i++)
    {
      wdomain[i]    = new double[2];
      wdomain[i][0] = w_min[w_indices[i]%3];
      wdomain[i][1] = w_max[w_indices[i]%3];
      VLOG(1) << "setting " << w_min[w_indices[i]%3] << " " << w_max[w_indices[i]%3];
    }
  }

  W->setDomains(wdomain);
  W->setBinSizes(wbins);

  VLOG(100) << "W:";
  VLOG(100) << *W;

  ULContainer* Wd = W->discretise();

  VLOG(100) << "Wd: " << endl << *Wd;

  // A data
  DContainer* a_domains = NULL;
  DContainer* A         = NULL;
  if(FLAGS_ai == "")
  {
    a_domains = csv->read(a_domain_file);
    A         = csv->read(a_states);
  }
  else
  {
    vector<int> a_indices = int_tokenizer(FLAGS_ai);
    a_domains = csv->read(a_domain_file, a_indices);
    A         = csv->read(a_states, a_indices);
  }
  A         = A->drop(1); // header

  check_domains("A", a_domains);
  VLOG(1) << "A domains:";
  VLOG(1) << *a_domains;
  VLOG(1) << "A states:";
  VLOG(1) << *A;

  double **adomain = new double*[A->columns()];
  int *abins = new int[A->columns()];
  for(int i = 0; i < A->columns(); i++)
  {
    adomain[i]    = new double[2];
    adomain[i][0] = (*a_domains)(0,i);
    adomain[i][1] = (*a_domains)(1,i);
    abins[i]      = FLAGS_abins;
  }

  A->setDomains(adomain);
  A->setBinSizes(abins);

  VLOG(100) << "A:";
  VLOG(100) << *A;

  ULContainer* Ad = A->discretise();

  ULContainer* W1 =  Wd->drop(-1);
  ULContainer* A1 =  Ad->drop(-1);
  ULContainer* W2 =  Wd->drop(1);


  if(FLAGS_csv)
  {
    ULContainer *Wdc = W->discretiseByColumn();
    csv->write(FLAGS_d + "/W.csv",   W);
    csv->write(FLAGS_d + "/Wd.csv",  Wd);
    csv->write(FLAGS_d + "/Wdc.csv", Wdc);
    csv->write(FLAGS_d + "/W1.csv",  W1);
    csv->write(FLAGS_d + "/W2.csv",  W2);
    csv->write(FLAGS_d + "/W_domains.csv", w_domains);
    ofstream o(FLAGS_d + "/W_min_max.csv");
    o << w_min[0] << "," << w_min[1] << "," << w_min[2] << endl;
    o << w_max[0] << "," << w_max[1] << "," << w_max[2] << endl;
    o.close();

    ULContainer *Adc = A->discretiseByColumn();
    csv->write(FLAGS_d + "/A.csv",   A);
    csv->write(FLAGS_d + "/A1.csv",  A1);
    csv->write(FLAGS_d + "/Ad.csv",  Ad);
    csv->write(FLAGS_d + "/Adc.csv", Adc);
    csv->write(FLAGS_d + "/A_domains.csv", a_domains);
  }
  

  double     mc  = entropy::sparse::MC_W(W2, W1, A1);
  DContainer* mcd = entropy::sparse::state::MC_W(W2, W1, A1);

  cout << "MC: " << mc << endl;

  stringstream sst;
  sst.str("");
  sst << FLAGS_d << "/mc_w-averaged.txt";
  VLOG(10) << "writing \"" << sst.str() << "\"";
  ofstream output;
#ifdef __APPLE__
  output.open(sst.str());
#else
  output.open(sst.str().c_str());
#endif
  output << mc << endl;
  output.close();

  sst.str("");
  sst << FLAGS_d << "/mc_w-state_dependent.csv";
  csv->write(sst.str().c_str(), mcd);

  VLOG(1) << "done.";
}
