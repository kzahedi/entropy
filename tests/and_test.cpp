#define BOOST_TEST_MODULE gis_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <glog/logging.h>

#define ITERATIONS  500000
#define DEBUG_LEVEL 0

#include <fstream>
#include <iostream>
#include <string>

#include <math.h>

# define EPSILON         0.01
# define ERROR_THRESHOLD 0.0001

#include <entropy++/Container.h>
#include <entropy++/iterativescaling/GIS.h>
#include <entropy++/iterativescaling/SCGIS.h>
#include <entropy++/iterativescaling/KL.h>


using namespace std;
using namespace entropy;
using namespace entropy::iterativescaling;

BOOST_AUTO_TEST_SUITE(AND)

BOOST_AUTO_TEST_CASE(AND_GIS)
{
  google::InitGoogleLogging("");
  FLAGS_logtostderr = 1;
  FLAGS_v = DEBUG_LEVEL;
  VLOG(10) << "Setting up X data";
  ULContainer *xData = new ULContainer(4,2);
  *xData << 0 << 0;
  *xData << 0 << 1;
  *xData << 1 << 0;
  *xData << 1 << 1;

  VLOG(10) << "Setting up Y data";
  ULContainer *yData = new ULContainer(4,1);
  *yData << 0;
  *yData << 0;
  *yData << 0;
  *yData << 1;

  ////////////////////////////////////////////////////////////////////////////////
  // independent model
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<int> > ia;
  vector<vector<int> > ib;

  vector<int> iaa;
  iaa.push_back(0);
  ia.push_back(iaa);
  vector<int> iab;
  iab.push_back(1);
  ia.push_back(iab);

  vector<int> ibb;
  ibb.push_back(0);
  ib.push_back(ibb);

  vector<Feature*> features;
  features.push_back(new Feature(0,0));
  features.push_back(new Feature(1,0));

  VLOG(10) << "Setting up independent model";
  GIS* independentModel = new GIS();
  independentModel->setData(xData, yData);
  independentModel->setFeatures(ia,ib,features);
  independentModel->init();

  for(int i = 0; i < ITERATIONS; i++)
  {
    independentModel->iterate();
    if(i % 100 == 0)
    {
      VLOG(1000) << *independentModel;
      VLOG(100)  << i << ": " << independentModel->error() << endl;
    }
    if(independentModel->error() < ERROR_THRESHOLD)
    {
      cout << "Stopping after " << i << " iterations." << endl;
      break;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // dependent model
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<int> > da;
  vector<vector<int> > db;

  vector<int> daa;
  daa.push_back(0);
  daa.push_back(1);
  da.push_back(daa);

  vector<int> dbb;
  dbb.push_back(0);
  db.push_back(dbb);

  vector<Feature*> dfeatures;
  dfeatures.push_back(new Feature(0,0));

  GIS* dependentModel = new GIS();
  dependentModel->setData(xData, yData);
  dependentModel->setFeatures(da,db,dfeatures);
  dependentModel->init();

  for(int i = 0; i < ITERATIONS; i++)
  {
    dependentModel->iterate();
    // if(i % 100 == 0) cout << i << ": " << dependentModel->error() << endl;
    if(dependentModel->error() < ERROR_THRESHOLD)
    {
      cout << "Stopping after " << i << " iterations." << endl;
      break;
    }

  }

  dependentModel->calculateProbabilities();

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
  cout << "AND-GIS (bits): " << kl->divergence2() << endl;
}

/*
BOOST_AUTO_TEST_CASE(AND_SCGIS)
{
  FLAGS_logtostderr = 1;
  FLAGS_v = DEBUG_LEVEL;
  VLOG(10) << "Setting up X data";
  ULContainer *xData = new ULContainer(4,2);
  *xData << 0 << 0;
  *xData << 0 << 1;
  *xData << 1 << 0;
  *xData << 1 << 1;

  VLOG(10) << "Setting up Y data";
  ULContainer *yData = new ULContainer(4,1);
  *yData << 0;
  *yData << 0;
  *yData << 0;
  *yData << 1;

  ////////////////////////////////////////////////////////////////////////////////
  // independent model
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<int> > ia;
  vector<vector<int> > ib;

  vector<int> iaa;
  iaa.push_back(0);
  ia.push_back(iaa);
  vector<int> iab;
  iab.push_back(1);
  ia.push_back(iab);

  vector<int> ibb;
  ibb.push_back(0);
  ib.push_back(ibb);

  vector<Feature*> features;
  features.push_back(new Feature(0,0));
  features.push_back(new Feature(1,0));

  VLOG(10) << "Setting up independent model";
  SCGIS* independentModel = new SCGIS();
  independentModel->setData(xData, yData);
  independentModel->setFeatures(ia,ib,features);
  independentModel->init();

  for(int i = 0; i < ITERATIONS; i++)
  {
    independentModel->iterate();
    // if(i % 100 == 0) cout << i << ": " << independentModel->error() << endl;
    VLOG(1000) << *independentModel;
    if(independentModel->error() < ERROR_THRESHOLD)
    {
      cout << "Stopping after " << i << " iterations." << endl;
      break;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // dependent model
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<int> > da;
  vector<vector<int> > db;

  vector<int> daa;
  daa.push_back(0);
  daa.push_back(1);
  da.push_back(daa);

  vector<int> dbb;
  dbb.push_back(0);
  db.push_back(dbb);

  vector<Feature*> dfeatures;
  dfeatures.push_back(new Feature(0,0));

  SCGIS* dependentModel = new SCGIS();
  dependentModel->setData(xData, yData);
  dependentModel->setFeatures(da,db,dfeatures);
  dependentModel->init();

  for(int i = 0; i < ITERATIONS; i++)
  {
    dependentModel->iterate();
    // if(i % 100 == 0) cout << i << ": " << dependentModel->error() << endl;
    if(dependentModel->error() < ERROR_THRESHOLD)
    {
      cout << "Stopping after " << i << " iterations." << endl;
      break;
    }

  }

  dependentModel->calculateProbabilities();

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
  cout << "AND-SCGIS (bits): " << kl->divergence2() << endl;
}
*/
BOOST_AUTO_TEST_SUITE_END()
