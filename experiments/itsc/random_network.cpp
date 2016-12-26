#include <iostream>
#include <fstream>

#include <stdlib.h>
#include <math.h>

#include <entropy++/Container.h>
#include <entropy++/Matrix.h>
#include <entropy++/Random.h>
#include <entropy++/iterativescaling/GIS.h>
#include <entropy++/iterativescaling/SCGIS.h>
#include <entropy++/iterativescaling/KL.h>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <string>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>


using namespace std;
using namespace boost;
using namespace entropy;
using namespace entropy::iterativescaling;

DEFINE_int64(i,  1000,  "nr. of iterations for each experiment");
DEFINE_int64(l,  100,   "nr. of last iterations to keep for the analysis");
DEFINE_int64(e,  100,   "nr. of experiments with random weight matrix");
DEFINE_int64(c,  100,   "nr. of convergence iterations");
DEFINE_double(t, 0.1,   "convergence abort threshold");
DEFINE_int64(n,  5,     "nr. of neurones");
DEFINE_double(b, 0.1,   "beta");
DEFINE_bool(g,   false, "use GIS instead of SCGIS");
DEFINE_string(o, "results.csv", "output file");

# define SIGM(x) (1.0 / (1.0 + exp(-x)))

void updateNetwork(Matrix& X, Matrix& W, double beta)
{
  Matrix Y(FLAGS_n, 1, 0.0);

  for(int i = 0; i < FLAGS_n; i++)
  {
    for(int j = 0; j < FLAGS_n; j++)
    {
      Y(i,0) += W(i,j) * X(j,0);
    }
  }

  Y *= (-2.0 * beta);

  for(int i = 0; i < FLAGS_n; i++) Y(i,0) = SIGM(Y(i,0));
  for(int i = 0; i < FLAGS_n; i++) Y(i,0) = (Random::unit() <= Y(i,0))?1:0;
  for(int i = 0; i < FLAGS_n; i++) X(i,0) = Y(i,0);

}

int main(int argc, char** argv)
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  FLAGS_log_dir = ".";
  FLAGS_logtostderr = true;
  Random::initialise();

  VLOG(0) << "Starting an experiment with: ";
  VLOG(0) << "  " << FLAGS_i << " iterations.";
  VLOG(0) << "  " << FLAGS_l << " last iterations will be used for the analysis.";
  VLOG(0) << "  " << FLAGS_e << " random initial conditions.";
  VLOG(0) << "  " << FLAGS_b << " temperature.";

  vector<double> results;

  for(int e = 0; e < FLAGS_e; e++)
  {
    Matrix X(FLAGS_n, 1);
    Matrix W(FLAGS_n, FLAGS_n);

    ULContainer data(FLAGS_i+1, FLAGS_n); 
    ULContainer first(1, FLAGS_n); 

    for(int i = 0; i < FLAGS_n; i++)
    {
      X(i,0) = (Random::unit() < 0.5)?0:1;
      data  << X(i,0);
      first << X(i,0);
      for(int j = 0; j < FLAGS_n; j++)
      {
        W(i,j) = Random::unit();
      }
    }

    for(int i = 0; i < FLAGS_i; i++)
    {
      updateNetwork(X, W, FLAGS_b);
      for(int j = 0; j < FLAGS_n; j++)
      {
        data << X(j,0);
      }
    }

    VLOG(10) << "Initial X:";
    VLOG(10) << X;
    VLOG(10) << "Weight Matrix:";
    VLOG(10) << W;

    VLOG(20) << "Original Data";
    VLOG(20) << data;

    ULContainer* input = data.drop(FLAGS_i - FLAGS_l);

    VLOG(10) << "Data for Analysis";
    VLOG(10) << *input;

    // cout << data << endl;

    vector<vector<int> > px;
    vector<vector<int> > py;

    vector<int> pa;
    vector<int> pb;
    for(int i = 0; i < FLAGS_n; i++)
    {
      pa.push_back(i);
      pb.push_back(i);
    }
    px.push_back(pa);
    py.push_back(pb);

    vector<vector<int> > qx;
    vector<vector<int> > qy;
    for(int i = 0; i < FLAGS_n; i++)
    {
      vector<int> qa;
      qa.push_back(i);
      qx.push_back(qa);

      vector<int> qb;
      qb.push_back(i);
      qy.push_back(qb);
    }

    ULContainer *Xt1 = input->drop(-1);
    ULContainer *Xt2 = input->drop(1);

    vector<Feature*> pfeatures;
    pfeatures.push_back(new Feature(0,0));

    vector<Feature*> qfeatures;
    for(int i = 0; i < FLAGS_n; i++)
    {
      qfeatures.push_back(new Feature(i, i));
    }

    IS* independentModel = NULL;
    if(FLAGS_g)
    {
      independentModel = new GIS();
    }
    else
    {
      independentModel = new SCGIS();
    }
    independentModel->setData(Xt1, Xt2);
    independentModel->setFeatures(qx,qy,qfeatures);
    independentModel->init();

    IS* dependentModel = NULL;
    if(FLAGS_g)
    {
      cout << "using GIS" << endl;
      dependentModel = new GIS();
    }
    else
    {
      cout << "using SCGIS" << endl;
      dependentModel = new SCGIS();
    }
    dependentModel->setData(Xt1, Xt2);
    dependentModel->setFeatures(px,py,pfeatures);
    dependentModel->init();

    double error = 0.0;
    int dependent_iterations = 0;
    int independent_iterations = 0;
    for(int i = 0; i < FLAGS_c; i++)
    {
      dependentModel->iterate();
      VLOG(11) << "Dependent model error after " << i << " iterations: " << dependentModel->error();
      if(dependentModel->error() < FLAGS_t)
      {
        dependent_iterations = i;
        break;
      }
    }
    for(int i = 0; i < FLAGS_c; i++)
    {
      independentModel->iterate();
      VLOG(11) << "Independent model error after " << i << " iterations: " << independentModel->error();
      if(independentModel->error() < FLAGS_t)
      {
        independent_iterations = i;
        break;
      }
    }

    KL* kl = new KL(dependentModel, independentModel);
    double r = kl->divergence();
    VLOG(0) << "Result: " << r;

    delete kl;
    if(FLAGS_g)
    {
      delete (GIS*)dependentModel;
      delete (GIS*)independentModel;
    }
    else
    {
      delete (SCGIS*)dependentModel;
      delete (SCGIS*)independentModel;
    }

    results.push_back(r);
  }

  stringstream sst;
  sst << "Results:";
  for(int i = 0; i < results.size(); i++) sst << " " << results[i];
  VLOG(0) << sst.str();

  ofstream out(FLAGS_o.c_str());
  out << FLAGS_b << ",";
  for(int e = 0; e < FLAGS_e - 1; e++)
  {
    out << results[e] << ",";
  }
  out << results[results.size() - 1];
  out << endl;
  out.close();
  VLOG(0) << "done.";
}
