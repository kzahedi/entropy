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
  vector<vector<int> > ib;

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
  for(int i = 0; i < 200; i++)
  {
    independentModel->iterate();
  }
  cout << (*independentModel) << endl;
  independentModel->calculateProbabilities();

  ////////////////////////////////////////////////////////////////////////////////
  // dependent model
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<int> > da;
  vector<vector<int> > db;

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
  for(int i = 0; i < 200; i++)
  {
    dependentModel->iterate();
  }
  cout << (*dependentModel) << endl;

  dependentModel->calculateProbabilities();
/*
  Matrix ipycx(2,4);
  ipycx(0,0) = 1.0;
  ipycx(0,1) = 1.0;
  ipycx(0,2) = 1.0;
  ipycx(1,3) = 1.0;

  Matrix ipx(1,4);
  ipx(0,0) = 1.0/4.0;
  ipx(0,1) = 1.0/4.0;
  ipx(0,2) = 1.0/4.0;
  ipx(0,3) = 1.0/4.0;


  BOOST_CHECK((fabs(ipycx(0,0) - dependentModel->p_y_c_x(0,0)) < EPSILON));
  BOOST_CHECK((fabs(ipycx(0,1) - dependentModel->p_y_c_x(0,1)) < EPSILON));
  BOOST_CHECK((fabs(ipycx(0,2) - dependentModel->p_y_c_x(0,2)) < EPSILON));
  BOOST_CHECK((fabs(ipycx(0,3) - dependentModel->p_y_c_x(0,3)) < EPSILON));
  BOOST_CHECK((fabs(ipycx(1,0) - dependentModel->p_y_c_x(1,0)) < EPSILON));
  BOOST_CHECK((fabs(ipycx(1,1) - dependentModel->p_y_c_x(1,1)) < EPSILON));
  BOOST_CHECK((fabs(ipycx(1,2) - dependentModel->p_y_c_x(1,2)) < EPSILON));
  BOOST_CHECK((fabs(ipycx(1,3) - dependentModel->p_y_c_x(1,3)) < EPSILON));
  BOOST_CHECK((fabs(ipx(0,0)   - dependentModel->p_x(0))       < EPSILON));
  BOOST_CHECK((fabs(ipx(0,1)   - dependentModel->p_x(1))       < EPSILON));
  BOOST_CHECK((fabs(ipx(0,2)   - dependentModel->p_x(2))       < EPSILON));
  BOOST_CHECK((fabs(ipx(0,3)   - dependentModel->p_x(3))       < EPSILON));

  ////////////////////////////////////////////////////////////////////////////////
  // Synergy
  ////////////////////////////////////////////////////////////////////////////////

  KL* kl = new KL(dependentModel, independentModel);
  cout << "AND (bits): " << kl->divergence2() << endl;
  cout << "AND (nats): " << kl->divergenceN() << endl; */
} }
#endif
