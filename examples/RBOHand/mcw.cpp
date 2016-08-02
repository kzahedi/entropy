#include "lib/Experiment.h"

#include <iostream>
#include <fstream>

#include <stdlib.h>

#include <entropy++/Csv.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/sparse/state/MC_W.h>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <string>

#include <btkAcquisitionFileReader.h>
#include <btkAcquisitionFileWriter.h>
#include <btkMacro.h> // btkStripPathMacro

#include <entropy++/Csv.h>

#define SET_DOMAIN_AND_BINS(a,d,b) \
  a->setDomains(d);   \
  a->setBinSizes(b);

using namespace std;

double x_min;
double x_max;

double y_min;
double y_max;

double z_min;
double z_max;

vector<double> amin;
vector<double> amax;

vector<string> files;

DEFINE_string(W,    "",   "W file");
DEFINE_string(A,    "",   "A file");
DEFINE_string(Wd,   "",   "W domain file");
DEFINE_string(Ad,   "",   "A domain file");
DEFINE_int64(wbins, 1000, "world state bins");
DEFINE_int64(abins, 100,  "action bins");

int main(int argc, char** argv)
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;

  VLOG(10) << "A File: " << FLAGS_W;
  VLOG(10) << "W File: " << FLAGS_A;
  VLOG(10) << "W Bins: " << FLAGS_wbins;
  VLOG(10) << "A Bins: " << FLAGS_abins;

  Csv*        csv = new Csv();
  Experiment* e   = new Experiment();

  e->setWorldStateBins(FLAGS_wbins);
  e->setActionStateBins(FLAGS_abins);

  e->wFile("", FLAGS_W);
  e->aFile("", FLAGS_A);

  double** w_domains = new double[e->wsize()];
  for(int i = 0; i < e->wsize(); i++)
  {
    w_domains[i] = new double[2];
    w_domains[i][0] = e->wmin(0,i);
    w_domains[i][1] = e->wmax(0,i);
  }




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
