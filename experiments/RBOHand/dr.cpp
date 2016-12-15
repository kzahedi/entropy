#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <math.h>

#include <entropy++/Csv.h>
#include <entropy++/Container.h>
#include <entropy++/MI.h>

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
using namespace entropy;

DEFINE_string(i, "",    "input file");
DEFINE_string(o, "",    "output file prefix");
DEFINE_double(m, 500.0, "max value for grasp distance");
DEFINE_int64(b,  30,    "nr. of bins");


int main(int argc, char** argv)
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  FLAGS_log_dir = ".";

  // data
  VLOG(10) << "reading " << FLAGS_i;
  DContainer* data = Csv::read(FLAGS_i);
  int columns = data->columns();
  VLOG(10) << "Data's dimensions: " << data->rows() << ", " << data->columns();

  VLOG(10) << "splitting the data ...";
  DContainer *mc    = data->columns(1,columns-2);
  DContainer *grasp = data->columns(1,columns-1);

  vector<int> indices;
  for(int i = 0; i < 93*93; i++) indices.push_back(i);
  DContainer *ndata = data->columns(indices);
  VLOG(10) << "Data's new dimensions: " << data->rows() << ", " << data->columns();
  delete data;

  ndata->setDomains(ndata->min(), ndata->max());
  ndata->setBinSizes(FLAGS_b);

  ULContainer* dd = ndata->discretiseByColumn();
  delete ndata;

  double mc_min = mc->min();
  double mc_max = mc->max();

  double grasp_min = grasp->min();
  double grasp_max = grasp->max();

  mc->setDomains(mc_min, mc_max);
  mc->setBinSizes(FLAGS_b);

  grasp->setDomains(grasp_min, grasp_max);
  grasp->setBinSizes(FLAGS_b);

  ULContainer* dmc = mc->discretise();
  delete mc;

  ULContainer* dgrasp = grasp->discretise();
  delete grasp;

  DContainer* r_mc    = new DContainer(93,93);
  DContainer* r_grasp = new DContainer(93,93);

  for(int i = 0; i < 93*93; i++)
  {
    ULContainer *d = dd->columns(1,i);
    (*r_mc) << entropy::MI(d, dmc);
    (*r_grasp) << entropy::MI(d, dgrasp);
  }

  string out_mc    = FLAGS_o + string(".mi_mc_w.csv");
  string out_grasp = FLAGS_o + string(".mi_grasp.csv");
  Csv::write(out_mc,    r_mc);
  Csv::write(out_grasp, r_grasp);
}

