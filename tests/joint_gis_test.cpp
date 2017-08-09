#define BOOST_TEST_MODULE joint_gis_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// #include <glog/logging.h>
// #include <boost/test/included/unit_test.hpp>

#define DEBUG_LEVEL 0
#define ITERATIONS 500000
#define TEST_AND

#include <entropy++/Container.h>
#include <entropy++/Csv.h>
#include <entropy++/iterativescaling/Model.h>
#include <entropy++/iterativescaling/GIS.h>
#include <entropy++/iterativescaling/Joint_GIS.h>
#include <entropy++/iterativescaling/KL.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/SparseMatrix.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <string>

#include <math.h>

# define EPSILON         0.01
# define ERROR_THRESHOLD 0.0001

# define TOLERANCE(a) boost::test_tools::tolerance(a)

# define BINS 30
// # define BINS 3


#ifndef PARENT
#define PARENT "/Users/zahedi/projects/entropy/experiments/hopping/data/"
#endif

# define MAX3(a,b,c) ENTROPY_MAX(a,ENTROPY_MAX(b,c))
# define MIN3(a,b,c) ENTROPY_MIN(a,ENTROPY_MIN(b,c))

# define POS 0
# define VEL 1
# define ACC 2
# define ACT 3
# define MSI 4 // muscle sensor input

# define POS_MATLAB_INDEX                 (2  - 1)
# define VEL_MATLAB_INDEX                 (3  - 1)
# define ACC_MATLAB_INDEX                 (4  - 1)
# define ACT_MATLAB_INDEX                 (10 - 1)
# define MUSCLE_SENSOR_INPUT_MATLAB_INDEX (5  - 1)



using namespace std;
using namespace entropy;
using namespace entropy::iterativescaling;

entropy::SparseMatrix pz;
entropy::SparseMatrix pyz;
entropy::SparseMatrix px_c_yz;
entropy::SparseMatrix pxyz;

BOOST_AUTO_TEST_SUITE(Logic)

#ifdef TEST_AND
BOOST_AUTO_TEST_CASE(AND)
{
 // google::InitGoogleLogging("");
 // FLAGS_logtostderr = 1;
 // FLAGS_v = DEBUG_LEVEL;
 // VLOG(10) << "Setting up X data";
  ULContainer *xData = new ULContainer(4,3);
  *xData << 0 << 0 << 0;
  *xData << 0 << 1 << 0;
  *xData << 1 << 0 << 0;
  *xData << 1 << 1 << 1;

  ULContainer *yData = new ULContainer(0,0);

  ULContainer *XAlphabet = new ULContainer(8,3);
  *XAlphabet << 0 << 0 << 0;
  *XAlphabet << 0 << 0 << 1;
  *XAlphabet << 0 << 1 << 0;
  *XAlphabet << 0 << 1 << 1;
  *XAlphabet << 1 << 0 << 0;
  *XAlphabet << 1 << 0 << 1;
  *XAlphabet << 1 << 1 << 0;
  *XAlphabet << 1 << 1 << 1;

  ////////////////////////////////////////////////////////////////////////////////
  // independent model
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<int> > ia;

  vector<int> iaa;
  iaa.push_back(0);
  iaa.push_back(2);
  ia.push_back(iaa);
  vector<int> iab;
  iab.push_back(1);
  iab.push_back(2);
  ia.push_back(iab);


  vector<Feature*> features;
  features.push_back(new Feature(0,-1));
  features.push_back(new Feature(1,-1));

//  VLOG(10) << "Setting up independent model";
  Joint_GIS* independentModel = new Joint_GIS();
  independentModel->setData(xData, yData);
  independentModel->setAlphabet(XAlphabet);
  independentModel->setFeatures(ia,features);
  independentModel->init();
  cout << "independent: " << endl;
  for(int i = 0; i < 2000; i++)
  {
    independentModel->iterate();
  }
  cout << (*independentModel) << endl;

  vector<double> q = independentModel->getProbabilities();

  ////////////////////////////////////////////////////////////////////////////////
  // dependent model
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<int> > da;

  vector<int> daa;
  daa.push_back(0);
  daa.push_back(1);
  daa.push_back(2);
  da.push_back(daa);

  vector<Feature*> dfeatures;
  dfeatures.push_back(new Feature(0,-1));

  Joint_GIS* dependentModel = new Joint_GIS();
  dependentModel->setData(xData, yData);
  dependentModel->setAlphabet(XAlphabet);
  dependentModel->setFeatures(da,dfeatures);
  dependentModel->init();

  cout << "dependent" << endl;
  for(int i = 0; i < 2000; i++)
  {
    dependentModel->iterate();
  }
  cout << (*dependentModel) << endl;

  cout << "KL: " <<  dependentModel->calculateKL(q) << endl;
}
}
BOOST_AUTO_TEST_CASE(ANDCMI)
{
  ULContainer *xData = new ULContainer(4,3);
  *xData << 0 << 0 << 0;
  *xData << 0 << 1 << 0;
  *xData << 1 << 0 << 0;
  *xData << 1 << 1 << 1;

  ULContainer *XAlphabet = new ULContainer(8,3);
  *XAlphabet << 0 << 0 << 0;
  *XAlphabet << 0 << 0 << 1;
  *XAlphabet << 0 << 1 << 0;
  *XAlphabet << 0 << 1 << 1;
  *XAlphabet << 1 << 0 << 0;
  *XAlphabet << 1 << 0 << 1;
  *XAlphabet << 1 << 1 << 0;
  *XAlphabet << 1 << 1 << 1;

  ULContainer *yData = new ULContainer(0,0);

  ////////////////////////////////////////////////////////////////////////////////
  // independent model
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<int> > ia;


  vector<int> iaa;
  iaa.push_back(0);
  iaa.push_back(2);
  ia.push_back(iaa);

  vector<Feature*> features;
  features.push_back(new Feature(0,-1));

  Joint_GIS* independentModel = new Joint_GIS();
  independentModel->setData(xData, yData);
  independentModel->setAlphabet(XAlphabet);
  independentModel->setFeatures(ia,features);
  independentModel->init();

  for(int i = 0; i < 500; i++)
  {
    independentModel->iterate();
  }
  cout << (*independentModel) << endl;
  vector<double> q = independentModel->getProbabilities();

  ////////////////////////////////////////////////////////////////////////////////
  // dependent model
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<int> > da;

  vector<int> daa;
  daa.push_back(0);
  daa.push_back(1);
  daa.push_back(2);
  da.push_back(daa);

  vector<Feature*> dfeatures;
  dfeatures.push_back(new Feature(0,-1));

  Joint_GIS* dependentModel = new Joint_GIS();
  dependentModel->setData(xData, yData);
  dependentModel->setAlphabet(XAlphabet);
  dependentModel->setFeatures(da, dfeatures);
  dependentModel->init();

  for(int i = 0; i < 500; i++)
  {
    dependentModel->iterate();
    if(dependentModel->error() < 0.000000001) break;
  }
  cout << (*dependentModel) << endl;

  cout << "KL: " <<  dependentModel->calculateKL(q) << endl;

}

#endif
