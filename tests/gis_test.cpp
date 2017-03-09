#include "gis_test.h"

#include <entropy++/Container.h>
#include <entropy++/Csv.h>
#include <entropy++/iterativescaling/Model.h>
#include <entropy++/iterativescaling/GIS.h>
#include <entropy++/iterativescaling/KL.h>
#include <entropy++/sparse/MC_W.h>

#include <fstream>
#include <iostream>
#include <string>

#include <math.h>

# define EPSILON         0.01
# define ERROR_THRESHOLD 0.01

# define BINS 3


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


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( gisTest );


void gisTest::testAND()
{
  ULContainer *xData = new ULContainer(8,3);
  *xData << 0 << 0 << 0;
  *xData << 0 << 1 << 0;
  *xData << 1 << 0 << 2;
  *xData << 1 << 1 << 0;
  *xData << 0 << 0 << 0;
  *xData << 0 << 1 << 0;
  *xData << 1 << 0 << 0;
  *xData << 1 << 1 << 1;

  ULContainer *yData = new ULContainer(8,2);
  *yData << 0 << 0;
  *yData << 0 << 0;
  *yData << 0 << 0;
  *yData << 1 << 0;
  *yData << 0 << 1;
  *yData << 0 << 0;
  *yData << 0 << 0;
  *yData << 1 << 1;

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

  GIS* independentModel = new GIS();
  independentModel->setData(xData, yData);
  independentModel->setFeatures(ia,ib,features);
  independentModel->init();

  for(int i = 0; i < 10; i++)
  {
    independentModel->iterate();
    if(independentModel->error() < 0.000000001) break;
  }

  independentModel->calculateProbabilities();

   cout << "Final Error: " << independentModel->error() << endl;
   cout << "AND: Independent model: " << endl << *independentModel << endl;

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

  for(int i = 0; i < 10; i++)
  {
    dependentModel->iterate();
     cout << "Error: " << dependentModel->error() << endl;
     cout << *dependentModel << endl;
    if(dependentModel->error() < 0.000000001) break;
  }

  dependentModel->calculateProbabilities();

   cout << "Final Error: " << dependentModel->error() << endl;
   cout << "AND: Dependent model: " << endl << *dependentModel << endl;

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

  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,0), dependentModel->p_y_c_x(0,0), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,1), dependentModel->p_y_c_x(0,1), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,2), dependentModel->p_y_c_x(0,2), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,3), dependentModel->p_y_c_x(0,3), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,0), dependentModel->p_y_c_x(1,0), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,1), dependentModel->p_y_c_x(1,1), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,2), dependentModel->p_y_c_x(1,2), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,3), dependentModel->p_y_c_x(1,3), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,0), dependentModel->p_x(0), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,1), dependentModel->p_x(1), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,2), dependentModel->p_x(2), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,3), dependentModel->p_x(3), EPSILON);

  ////////////////////////////////////////////////////////////////////////////////
  // Synergy
  ////////////////////////////////////////////////////////////////////////////////

  KL* kl = new KL(dependentModel, independentModel);
  // CPPUNIT_ASSERT_DOUBLES_EQUAL(EPSILON, kl->divergence(), EPSILON);
  cout << "AND (bits): " << kl->divergence2() << endl;
  cout << "AND (nats): " << kl->divergenceN() << endl;
}

void gisTest::testANDCMI()
{
  ULContainer *xData = new ULContainer(8,3);
    *xData << 0 << 0 << 0;
    *xData << 0 << 1 << 0;
    *xData << 1 << 0 << 1;
    *xData << 1 << 1 << 0;
    *xData << 0 << 0 << 1;
    *xData << 0 << 1 << 1;
    *xData << 1 << 0 << 1;
    *xData << 1 << 1 << 1;

    ULContainer *yData = new ULContainer(8,2);
    *yData << 0 << 0;
    *yData << 0 << 0;
    *yData << 0 << 0;
    *yData << 1 << 0;
    *yData << 0 << 1;
    *yData << 0 << 1;
    *yData << 0 << 1;
    *yData << 1 << 1;

    ////////////////////////////////////////////////////////////////////////////////
    // independent model
    ////////////////////////////////////////////////////////////////////////////////
    vector<vector<int> > ia;
    vector<vector<int> > ib;

    vector<int> iaa;
    iaa.push_back(1);
    ia.push_back(iaa);

    vector<int> ibb;
    ibb.push_back(0);
    ib.push_back(ibb);

    vector<Feature*> features;
    features.push_back(new Feature(0,0));

    GIS* independentModel = new GIS();
    independentModel->setData(xData, yData);
    independentModel->setFeatures(ia,ib,features);
    independentModel->init();

    for(int i = 0; i < 11; i++)
    {
      independentModel->iterate();
      if(independentModel->error() < 0.000000001) break;
    }

    independentModel->calculateProbabilities();

    // cout << "Final Error: " << independentModel->error() << endl;
    // cout << "AND: Independent model: " << endl << *independentModel << endl;

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

    for(int i = 0; i < 11; i++)
    {
      dependentModel->iterate();
      // cout << "Error: " << dependentModel->error() << endl;
      // cout << *dependentModel << endl;
      if(dependentModel->error() < 0.000000001) break;
    }

    dependentModel->calculateProbabilities();

    // cout << "Final Error: " << dependentModel->error() << endl;
    // cout << "AND: Dependent model: " << endl << *dependentModel << endl;

    Matrix dpycx(2,4,0.0);
    dpycx(0,0) = 1.0;
    dpycx(0,1) = 1.0;
    dpycx(0,2) = 1.0;
    dpycx(1,3) = 1.0;

    Matrix dpx(1,4, 0.0);
    dpx(0,0) = 1.0/4.0;
    dpx(0,1) = 1.0/4.0;
    dpx(0,2) = 1.0/4.0;
    dpx(0,3) = 1.0/4.0;

    Matrix ipycx(2,2);
    ipycx(0,0) = 1.0;
    ipycx(1,0) = 0.0;
    ipycx(0,1) = 0.5;
    ipycx(1,1) = 0.5;

    Matrix ipx(1,2, 0.0);
    ipx(0,0) = 1.0/2.0;
    ipx(0,1) = 1.0/2.0;

    CPPUNIT_ASSERT_EQUAL(2, independentModel->p_y_c_x()->rows());
    CPPUNIT_ASSERT_EQUAL(2, independentModel->p_y_c_x()->cols());
    CPPUNIT_ASSERT_EQUAL(2, independentModel->p_x()->rows());
    CPPUNIT_ASSERT_EQUAL(1, independentModel->p_x()->cols());

    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpycx(0,0), dependentModel->p_y_c_x(0,0), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpycx(0,1), dependentModel->p_y_c_x(0,1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpycx(0,2), dependentModel->p_y_c_x(0,2), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpycx(0,3), dependentModel->p_y_c_x(0,3), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpycx(1,0), dependentModel->p_y_c_x(1,0), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpycx(1,1), dependentModel->p_y_c_x(1,1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpycx(1,2), dependentModel->p_y_c_x(1,2), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpycx(1,3), dependentModel->p_y_c_x(1,3), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpx(0,0),   dependentModel->p_x(0),       EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpx(0,1),   dependentModel->p_x(1),       EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpx(0,2),   dependentModel->p_x(2),       EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(dpx(0,3),   dependentModel->p_x(3),       EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,0), independentModel->p_y_c_x(0,0), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,1), independentModel->p_y_c_x(0,1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,0), independentModel->p_y_c_x(1,0), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,1), independentModel->p_y_c_x(1,1), EPSILON);
   // CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,0),   independentModel->p_x(0),       EPSILON);
   // CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,1),   independentModel->p_x(1),       EPSILON);

    ////////////////////////////////////////////////////////////////////////////////
    // Synergy
    ////////////////////////////////////////////////////////////////////////////////

    KL* kl = new KL(dependentModel, independentModel);
    // CPPUNIT_ASSERT_DOUBLES_EQUAL(EPSILON, kl->divergence(), EPSILON);
    cout << "AND CMI (bits): " << kl->divergence2() << endl;
    cout << "AND CMI (nats): " << kl->divergenceN() << endl;
}

void gisTest::testOR()
{
  ULContainer *xData = new ULContainer(4,2);
  *xData << 0 << 0;
  *xData << 0 << 1;
  *xData << 1 << 0;
  *xData << 1 << 1;

  ULContainer *yData = new ULContainer(4,1);
  *yData << 0;
  *yData << 1;
  *yData << 1;
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

  GIS* independentModel = new GIS();
  independentModel->setData(xData, yData);
  independentModel->setFeatures(ia,ib,features);
  independentModel->init();

  for(int i = 0; i < 10; i++)
  {
    independentModel->iterate();
    if(independentModel->error() < ERROR_THRESHOLD) break;
  }

  independentModel->calculateProbabilities();

   cout << "OR: Independent model: " << endl << *independentModel << endl;

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

  for(int i = 0; i < 10; i++)
  {
    dependentModel->iterate();
     cout << "Error: " << dependentModel->error() << endl;
     cout << *dependentModel << endl;
    if(dependentModel->error() < ERROR_THRESHOLD) break;
  }

  dependentModel->calculateProbabilities();

   cout << "Final Error: " << dependentModel->error() << endl;
   cout << "AND: Dependent model: " << endl << *dependentModel << endl;

  Matrix ipycx(2,4);
  ipycx(0,0) = 1.0;
  ipycx(1,1) = 1.0;
  ipycx(1,2) = 1.0;
  ipycx(1,3) = 1.0;

  Matrix ipx(1,4);
  ipx(0,0) = 1.0/4.0;
  ipx(0,1) = 1.0/4.0;
  ipx(0,2) = 1.0/4.0;
  ipx(0,3) = 1.0/4.0;

  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,0), dependentModel->p_y_c_x(0,0), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,1), dependentModel->p_y_c_x(0,1), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,2), dependentModel->p_y_c_x(0,2), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,3), dependentModel->p_y_c_x(0,3), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,0), dependentModel->p_y_c_x(1,0), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,1), dependentModel->p_y_c_x(1,1), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,2), dependentModel->p_y_c_x(1,2), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,3), dependentModel->p_y_c_x(1,3), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,0), dependentModel->p_x(0), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,1), dependentModel->p_x(1), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,2), dependentModel->p_x(2), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,3), dependentModel->p_x(3), EPSILON);


  ////////////////////////////////////////////////////////////////////////////////
  // Synergy
  ////////////////////////////////////////////////////////////////////////////////

  KL* kl = new KL(dependentModel, independentModel);
  // CPPUNIT_ASSERT_DOUBLES_EQUAL(EPSILON, kl->divergence(), EPSILON);
  cout << "OR (bits): " << kl->divergence2() << endl;
  cout << "OR (nats): " << kl->divergenceN() << endl;
}

void gisTest::testXOR()
{
  ULContainer *xData = new ULContainer(8,3);
  *xData << 0 << 0 << 0;
  *xData << 0 << 1 << 0;
  *xData << 1 << 0 << 0;
  *xData << 1 << 1 << 0;
  *xData << 0 << 0 << 1;
  *xData << 0 << 1 << 1;
  *xData << 1 << 0 << 1;
  *xData << 1 << 1 << 1;

  ULContainer *yData = new ULContainer(8,1);
  *yData << 0;
  *yData << 1;
  *yData << 1;
  *yData << 1;
  *yData << 1;
  *yData << 1;
  *yData << 1;
  *yData << 0;

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

  vector<int> iac;
  iac.push_back(2);
  ia.push_back(iac);

  vector<int> ibb;
  ibb.push_back(0);
  ib.push_back(ibb);

  vector<Feature*> features;
  features.push_back(new Feature(0,0));
  features.push_back(new Feature(1,0));
  features.push_back(new Feature(2,0));

  GIS* independentModel = new GIS();
  independentModel->setData(xData, yData);
  independentModel->setFeatures(ia,ib,features);
  independentModel->init();

  for(int i = 0; i < 10; i++)
  {
    independentModel->iterate();
    if(independentModel->error() < ERROR_THRESHOLD) break;
  }

  independentModel->calculateProbabilities();

   cout << "XOR: Independent model: " << endl << *independentModel << endl;

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

  vector<int> dbb;
  dbb.push_back(0);
  db.push_back(dbb);

  vector<Feature*> dfeatures;
  dfeatures.push_back(new Feature(0,0));

  GIS* dependentModel = new GIS();
  dependentModel->setData(xData, yData);
  dependentModel->setFeatures(da,db,dfeatures);
  dependentModel->init();

  for(int i = 0; i < 1000; i++)
  {
    dependentModel->iterate();
    independentModel->iterate();
     cout << "Error: " << dependentModel->error() << endl;
     cout << "Error: " << independentModel->error() << endl;
   //  cout << *dependentModel << endl;
   // if(dependentModel->error() < 0.000000001) break;
     KL* kl = new KL(dependentModel, independentModel);
     cout << " kl: " << kl->divergence2() << endl;
  }

  dependentModel->calculateProbabilities();

   cout << "Final Error: " << dependentModel->error() << endl;
   cout << "AND: Dependent model: " << endl << *dependentModel << endl;
/*
  Matrix ipycx(2,4);
  ipycx(0,0) = 1.0;
  ipycx(1,1) = 1.0;
  ipycx(1,2) = 1.0;
  ipycx(0,3) = 1.0;

  Matrix ipx(1,4);
  ipx(0,0) = 1.0/4.0;
  ipx(0,1) = 1.0/4.0;
  ipx(0,2) = 1.0/4.0;
  ipx(0,3) = 1.0/4.0;

  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,0), dependentModel->p_y_c_x(0,0), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,1), dependentModel->p_y_c_x(0,1), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,2), dependentModel->p_y_c_x(0,2), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(0,3), dependentModel->p_y_c_x(0,3), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,0), dependentModel->p_y_c_x(1,0), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,1), dependentModel->p_y_c_x(1,1), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,2), dependentModel->p_y_c_x(1,2), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipycx(1,3), dependentModel->p_y_c_x(1,3), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,0), dependentModel->p_x(0), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,1), dependentModel->p_x(1), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,2), dependentModel->p_x(2), EPSILON);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ipx(0,3), dependentModel->p_x(3), EPSILON);
*/

  ////////////////////////////////////////////////////////////////////////////////
  // Synergy
  ////////////////////////////////////////////////////////////////////////////////

  KL* kl = new KL(dependentModel, independentModel);
  // CPPUNIT_ASSERT_DOUBLES_EQUAL(EPSILON, kl->divergence(), EPSILON);
  cout << "XOR (bits): " << kl->divergence2() << endl;
  cout << "XOR (nats): " << kl->divergenceN() << endl;
}














void gisTest::testUnique()
{ /*
  cout << PARENT << "/dcmot.csv" << endl;

  Csv *csv = new Csv();
  DContainer *dcmot  = Csv::read(string(PARENT) + string("/dcmot.csv"),  4,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX);

  DContainer *muslin = Csv::read(string(PARENT) + string("/muslin.csv"), 5,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX,
                                 MUSCLE_SENSOR_INPUT_MATLAB_INDEX);

  DContainer *musfib = Csv::read(string(PARENT) + string("/musfib.csv"), 5,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX,
                                 MUSCLE_SENSOR_INPUT_MATLAB_INDEX);

  DContainer *dcmotW  = dcmot->columns(3, 0, 1, 2);
  DContainer *dcmotA  = dcmot->columns(1, 3);
  delete dcmot;

  DContainer *muslinW = muslin->columns(3, 0, 1, 2);
  DContainer *muslinA = muslin->columns(1, 3);
  delete muslin;

  DContainer *musfibW = musfib->columns(3, 0, 1, 2);
  DContainer *musfibA = musfib->columns(1, 3);
  delete musfib;

  double p_min = MIN3(dcmotW->min(0), muslinW->min(0), musfibW->min(0));
  double p_max = MAX3(dcmotW->max(0), muslinW->max(0), musfibW->max(0));

  double v_min = MIN3(dcmotW->min(1), muslinW->min(1), musfibW->min(1));
  double v_max = MAX3(dcmotW->max(1), muslinW->max(1), musfibW->max(1));

  double a_min = MIN3(dcmotW->min(2), muslinW->min(2), musfibW->min(2));
  double a_max = MAX3(dcmotW->max(2), muslinW->max(2), musfibW->max(2));

  double ac_min = dcmotA->min(0);
  double ac_max = dcmotA->max(0);

  cout << "Domains:" << endl;
  cout << "  Position:          " << p_min  << " " << p_max  << endl;
  cout << "  Velocity:          " << v_min  << " " << v_max  << endl;
  cout << "  Acceleration:      " << a_min  << " " << a_max  << endl;
  cout << "  Action (DC Motor): " << ac_min << " " << ac_max << endl;
  cout << "  Bins:              " << BINS   << endl;

  // normalising DCMot
  dcmotW->normaliseColumn(0, p_min,  p_max);
  dcmotW->normaliseColumn(1, v_min,  v_max);
  dcmotW->normaliseColumn(2, a_min,  a_max);

  dcmotA->normaliseColumn(0, ac_min, ac_max);

  // normalising MusLin
  muslinW->normaliseColumn(0, p_min,  p_max);
  muslinW->normaliseColumn(1, v_min,  v_max);
  muslinW->normaliseColumn(2, a_min,  a_max);

  // normalising MusFib
  musfibW->normaliseColumn(0, p_min,  p_max);
  musfibW->normaliseColumn(1, v_min,  v_max);
  musfibW->normaliseColumn(2, a_min,  a_max);

  double **wdomain = new double*[3];
  for(int i = 0; i < 3; i++)
  {
    wdomain[i] = new double[2];
    wdomain[i][0] = 0.0;
    wdomain[i][1] = 1.0;
  }
  int *wbins = new int[3];
  wbins[0]   = BINS;
  wbins[1]   = BINS;
  wbins[2]   = BINS;

  double **adomain = new double*[1];
  for(int i = 0; i < 1; i++)
  {
    adomain[i] = new double[2];
    adomain[i][0] = 0.0;
    adomain[i][1] = 1.0;
  }
  int *abins = new int[1];
  abins[0]   = BINS;

  dcmotW->setDomains(wdomain);
  dcmotW->setBinSizes(wbins);
  dcmotA->setDomains(adomain);
  dcmotA->setBinSizes(abins);

  muslinW->setDomains(wdomain);
  muslinW->setBinSizes(wbins);
  muslinA->setDomains(adomain);
  muslinA->setBinSizes(abins);

  musfibW->setDomains(wdomain);
  musfibW->setBinSizes(wbins);
  musfibA->setDomains(adomain);
  musfibA->setBinSizes(abins);

  cout << "discretising dcmot W" << endl;
  ULContainer *DdcmotW  = dcmotW->discretiseByColumn(); delete dcmotW;
  cout << "discretising dcmot A" << endl;
  ULContainer *DdcmotA  = dcmotA->discretiseByColumn(); delete dcmotA;
  ULContainer *dcW2 = DdcmotW->drop(1);
  ULContainer *dcW1 = DdcmotW->drop(-1);
  ULContainer *dcA1 = DdcmotA->drop(-1);

  cout << "discretising muslin W" << endl;
  ULContainer *DmuslinW  = muslinW->discretiseByColumn(); delete muslinW;
  cout << "discretising muslin A" << endl;
  ULContainer *DmuslinA  = muslinA->discretiseByColumn(); delete muslinA;
  ULContainer *mlW2 = DmuslinW->drop(1);
  ULContainer *mlW1 = DmuslinW->drop(-1);
  ULContainer *mlA1 = DmuslinA->drop(-1);

  cout << "discretising musfib W" << endl;
  ULContainer *DmusfibW  = musfibW->discretiseByColumn(); delete musfibW;
  cout << "discretising musfib A" << endl;
  ULContainer *DmusfibA  = musfibA->discretiseByColumn(); delete musfibA;
  ULContainer *mfW2 = DmusfibW->drop(1);
  ULContainer *mfW1 = DmusfibW->drop(-1);
  ULContainer *mfA1 = DmusfibA->drop(-1);


  vector<vector<int> > a;
  vector<vector<int> > b;

  vector<int> aa;
  aa.push_back(0);
  aa.push_back(1);
  aa.push_back(2);
  a.push_back(aa);

  vector<int> bb;
  bb.push_back(0);
  b.push_back(bb);

  vector<Feature*> features;
  features.push_back(new Feature(0,0));

  // a.Xindices.push_back(0);
  // a.Xindices.push_back(1);
  // a.Xindices.push_back(2);
  // a.Yindices.push_back(0);
  // mi.push_back(a);

  ////////////////////////////////////////////////////////////////////////////////
  // DC MOTOR
  ////////////////////////////////////////////////////////////////////////////////
  ULContainer *dcZ  = dcW2;
  ULContainer *dcXY = dcW1;

  GIS* dcgis = new GIS();
  dcgis->setData(dcXY, dcZ);
  dcgis->setFeatures(a,b,features);
  dcgis->init();

  CPPUNIT_ASSERT(dcgis->X()->rows() > 0);
  CPPUNIT_ASSERT(dcgis->X()->rows() == dcXY->rows());
  CPPUNIT_ASSERT(dcgis->Y()->rows() > 0);
  CPPUNIT_ASSERT(dcgis->Y()->rows() == dcZ->rows());

  // CPPUNIT_ASSERT(dcgis->uniqueX(0)->rows() < dcgis->X()->rows());
  // CPPUNIT_ASSERT(dcgis->uniqueY(0)->rows() < dcgis->Y()->rows());

  for(int i = 0; i < 500; i++)
  {
    dcgis->iterate();
    if(dcgis->error() < 0.00000001)
    {
      cout << "Error: " << dcgis->error() << " with " << i << " iterations" << endl;
      break;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////
  // MusFib MOTOR
  ////////////////////////////////////////////////////////////////////////////////
  ULContainer *mfZ  = mfW2;
  ULContainer *mfXY = mfW1;

  GIS* mfgis = new GIS();
  mfgis->setData(mfXY, mfZ);
  mfgis->setFeatures(a,b,features);
  mfgis->init();

  CPPUNIT_ASSERT(mfgis->X()->rows() > 0);
  CPPUNIT_ASSERT(mfgis->X()->rows() == mfXY->rows());
  CPPUNIT_ASSERT(mfgis->Y()->rows() > 0);
  CPPUNIT_ASSERT(mfgis->Y()->rows() == mfZ->rows());

  // CPPUNIT_ASSERT(mfgis->uniqueX(0)->rows() < mfgis->X()->rows());
  // CPPUNIT_ASSERT(mfgis->uniqueY(0)->rows() < mfgis->Y()->rows());

  for(int i = 0; i < 500; i++)
  {
    mfgis->iterate();
    if(mfgis->error() < 0.00000001)
    {
      cout << "Error: " << mfgis->error() << " with " << i << " iterations" << endl;
      break;
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // MusLin MOTOR
  ////////////////////////////////////////////////////////////////////////////////
  ULContainer *mlZ  = mlW2;
  ULContainer *mlXY = mlW1;

  GIS* mlgis = new GIS();
  mlgis->setData(mlXY, mlZ);
  mlgis->setFeatures(a,b,features);
  mlgis->init();

  CPPUNIT_ASSERT(mlgis->X()->rows() > 0);
  CPPUNIT_ASSERT(mlgis->X()->rows() == mlXY->rows());
  CPPUNIT_ASSERT(mlgis->Y()->rows() > 0);
  CPPUNIT_ASSERT(mlgis->Y()->rows() == mlZ->rows());

  // CPPUNIT_ASSERT(mlgis->uniqueX(0)->rows() < mlgis->X()->rows());
  // CPPUNIT_ASSERT(mlgis->uniqueY(0)->rows() < mlgis->Y()->rows());


  for(int i = 0; i < 500; i++)
  {
    mlgis->iterate();
    if(mlgis->error() < 0.00000001)
    {
      cout << "Error: " << mlgis->error() << " with " << i << " iterations" << endl;
      break;
    }
  } */
}

void gisTest::testUnique2()
{/*
  cout << PARENT << "/dcmot.csv" << endl;

  Csv *csv = new Csv();
  DContainer *dcmot  = Csv::read(string(PARENT) + string("/dcmot.csv"),  4,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX);

  DContainer *muslin = Csv::read(string(PARENT) + string("/muslin.csv"), 5,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX,
                                 MUSCLE_SENSOR_INPUT_MATLAB_INDEX);

  DContainer *musfib = Csv::read(string(PARENT) + string("/musfib.csv"), 5,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX,
                                 MUSCLE_SENSOR_INPUT_MATLAB_INDEX);

  DContainer *dcmotW  = dcmot->columns(3, 0, 1, 2);
  DContainer *dcmotA  = dcmot->columns(1, 3);
  delete dcmot;

  DContainer *muslinW = muslin->columns(3, 0, 1, 2);
  DContainer *muslinA = muslin->columns(1, 3);
  delete muslin;

  DContainer *musfibW = musfib->columns(3, 0, 1, 2);
  DContainer *musfibA = musfib->columns(1, 3);
  delete musfib;

  double p_min = MIN3(dcmotW->min(0), muslinW->min(0), musfibW->min(0));
  double p_max = MAX3(dcmotW->max(0), muslinW->max(0), musfibW->max(0));

  double v_min = MIN3(dcmotW->min(1), muslinW->min(1), musfibW->min(1));
  double v_max = MAX3(dcmotW->max(1), muslinW->max(1), musfibW->max(1));

  double a_min = MIN3(dcmotW->min(2), muslinW->min(2), musfibW->min(2));
  double a_max = MAX3(dcmotW->max(2), muslinW->max(2), musfibW->max(2));

  double ac_min = dcmotA->min(0);
  double ac_max = dcmotA->max(0);

  cout << "Domains:" << endl;
  cout << "  Position:          " << p_min  << " " << p_max  << endl;
  cout << "  Velocity:          " << v_min  << " " << v_max  << endl;
  cout << "  Acceleration:      " << a_min  << " " << a_max  << endl;
  cout << "  Action (DC Motor): " << ac_min << " " << ac_max << endl;
  cout << "  Bins:              " << BINS   << endl;

  // normalising DCMot
  dcmotW->normaliseColumn(0, p_min,  p_max);
  dcmotW->normaliseColumn(1, v_min,  v_max);
  dcmotW->normaliseColumn(2, a_min,  a_max);

  dcmotA->normaliseColumn(0, ac_min, ac_max);

  // normalising MusLin
  muslinW->normaliseColumn(0, p_min,  p_max);
  muslinW->normaliseColumn(1, v_min,  v_max);
  muslinW->normaliseColumn(2, a_min,  a_max);

  // normalising MusFib
  musfibW->normaliseColumn(0, p_min,  p_max);
  musfibW->normaliseColumn(1, v_min,  v_max);
  musfibW->normaliseColumn(2, a_min,  a_max);

  double **wdomain = new double*[3];
  for(int i = 0; i < 3; i++)
  {
    wdomain[i] = new double[2];
    wdomain[i][0] = 0.0;
    wdomain[i][1] = 1.0;
  }
  int *wbins = new int[3];
  wbins[0]   = BINS;
  wbins[1]   = BINS;
  wbins[2]   = BINS;

  double **adomain = new double*[1];
  for(int i = 0; i < 1; i++)
  {
    adomain[i] = new double[2];
    adomain[i][0] = 0.0;
    adomain[i][1] = 1.0;
  }
  int *abins = new int[1];
  abins[0]   = BINS;

  dcmotW->setDomains(wdomain);
  dcmotW->setBinSizes(wbins);
  dcmotA->setDomains(adomain);
  dcmotA->setBinSizes(abins);

  muslinW->setDomains(wdomain);
  muslinW->setBinSizes(wbins);
  muslinA->setDomains(adomain);
  muslinA->setBinSizes(abins);

  musfibW->setDomains(wdomain);
  musfibW->setBinSizes(wbins);
  musfibA->setDomains(adomain);
  musfibA->setBinSizes(abins);

  cout << "discretising dcmot W" << endl;
  ULContainer *DdcmotW  = dcmotW->discretiseByColumn(); delete dcmotW;
  cout << "discretising dcmot A" << endl;
  ULContainer *DdcmotA  = dcmotA->discretiseByColumn(); delete dcmotA;
  ULContainer *dcW2 = DdcmotW->drop(1);
  ULContainer *dcW1 = DdcmotW->drop(-1);
  ULContainer *dcA1 = DdcmotA->drop(-1);

  cout << "discretising muslin W" << endl;
  ULContainer *DmuslinW  = muslinW->discretiseByColumn(); delete muslinW;
  cout << "discretising muslin A" << endl;
  ULContainer *DmuslinA  = muslinA->discretiseByColumn(); delete muslinA;
  ULContainer *mlW2 = DmuslinW->drop(1);
  ULContainer *mlW1 = DmuslinW->drop(-1);
  ULContainer *mlA1 = DmuslinA->drop(-1);

  cout << "discretising musfib W" << endl;
  ULContainer *DmusfibW  = musfibW->discretiseByColumn(); delete musfibW;
  cout << "discretising musfib A" << endl;
  ULContainer *DmusfibA  = musfibA->discretiseByColumn(); delete musfibA;
  ULContainer *mfW2 = DmusfibW->drop(1);
  ULContainer *mfW1 = DmusfibW->drop(-1);
  ULContainer *mfA1 = DmusfibA->drop(-1);

  vector<vector<int> > a;
  vector<vector<int> > b;

  vector<int> aa;
  aa.push_back(0);
  vector<int> ab;
  ab.push_back(1);
  vector<int> ac;
  ac.push_back(2);
  a.push_back(aa);
  a.push_back(ab);
  a.push_back(ac);

  vector<int> bb;
  bb.push_back(0);
  b.push_back(bb);

  vector<Feature*> features;
  features.push_back(new Feature(0,0));
  features.push_back(new Feature(1,0));
  features.push_back(new Feature(2,0));

  ////////////////////////////////////////////////////////////////////////////////
  // DC MOTOR
  ////////////////////////////////////////////////////////////////////////////////
  ULContainer *dcZ  = dcW2;
  ULContainer *dcXY = dcW1;

  GIS* dcgis = new GIS();
  dcgis->setData(dcXY, dcZ);
  dcgis->setFeatures(a,b,features);
  dcgis->init();

  CPPUNIT_ASSERT(dcgis->X()->rows() > 0);
  CPPUNIT_ASSERT(dcgis->X()->rows() == dcXY->rows());
  CPPUNIT_ASSERT(dcgis->Y()->rows() > 0);
  CPPUNIT_ASSERT(dcgis->Y()->rows() == dcZ->rows());

  // CPPUNIT_ASSERT(dcgis->uniqueX(0)->rows() < dcgis->X()->rows());
  // CPPUNIT_ASSERT(dcgis->uniqueY(0)->rows() < dcgis->Y()->rows());

  for(int i = 0; i < 500; i++)
  {
    dcgis->iterate();
    if(dcgis->error() < 0.00000001)
    {
      cout << "Error: " << dcgis->error() << " with " << i << " iterations" << endl;
      break;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////
  // MusFib MOTOR
  ////////////////////////////////////////////////////////////////////////////////
  ULContainer *mfZ  = mfW2;
  ULContainer *mfXY = mfW1;

  GIS* mfgis = new GIS();
  mfgis->setData(mfXY, mfZ);
  mfgis->setFeatures(a,b,features);
  mfgis->init();

  CPPUNIT_ASSERT(mfgis->X()->rows() > 0);
  CPPUNIT_ASSERT(mfgis->X()->rows() == mfXY->rows());
  CPPUNIT_ASSERT(mfgis->Y()->rows() > 0);
  CPPUNIT_ASSERT(mfgis->Y()->rows() == mfZ->rows());

  // CPPUNIT_ASSERT(mfgis->uniqueX(0)->rows() < mfgis->X()->rows());
  // CPPUNIT_ASSERT(mfgis->uniqueY(0)->rows() < mfgis->Y()->rows());

  for(int i = 0; i < 500; i++)
  {
    mfgis->iterate();
    if(mfgis->error() < 0.00000001)
    {
      cout << "Error: " << mfgis->error() << " with " << i << " iterations" << endl;
      break;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////
  // MusLin MOTOR
  ////////////////////////////////////////////////////////////////////////////////
  ULContainer *mlZ  = mlW2;
  ULContainer *mlXY = mlW1;

  GIS* mlgis = new GIS();
  mlgis->setData(mlXY, mlZ);
  mlgis->setFeatures(a,b,features);
  mlgis->init();

  CPPUNIT_ASSERT(mlgis->X()->rows() > 0);
  CPPUNIT_ASSERT(mlgis->X()->rows() == mlXY->rows());
  CPPUNIT_ASSERT(mlgis->Y()->rows() > 0);
  CPPUNIT_ASSERT(mlgis->Y()->rows() == mlZ->rows());

  // CPPUNIT_ASSERT(mlgis->uniqueX(0)->rows() < mlgis->X()->rows());
  // CPPUNIT_ASSERT(mlgis->uniqueY(0)->rows() < mlgis->Y()->rows());

  for(int i = 0; i < 500; i++)
  {
    mlgis->iterate();
    if(mlgis->error() < 0.00000001)
    {
      cout << "Error: " << mlgis->error() << " with " << i << " iterations" << endl;
      break;
    }
  }
 */
}

void gisTest::testMC_W()
{
  cout << PARENT << "/dcmot_small.csv" << endl;

  Csv *csv = new Csv();
  DContainer *dcmot  = Csv::read(string(PARENT) + string("/dcmot.csv"),  4,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX);

  DContainer *muslin = Csv::read(string(PARENT) + string("/muslin.csv"), 5,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX,
                                 MUSCLE_SENSOR_INPUT_MATLAB_INDEX);

  DContainer *musfib = Csv::read(string(PARENT) + string("/musfib.csv"), 5,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX,
                                 MUSCLE_SENSOR_INPUT_MATLAB_INDEX);

  DContainer *dcmotW  = dcmot->columns(3, 0, 1, 2);
  DContainer *dcmotA  = dcmot->columns(1, 3);
  delete dcmot;

  DContainer *muslinW = muslin->columns(3, 0, 1, 2);
  DContainer *muslinA = muslin->columns(1, 3);
  delete muslin;

  DContainer *musfibW = musfib->columns(3, 0, 1, 2);
  DContainer *musfibA = musfib->columns(1, 3);
  delete musfib;

  double p_min = MIN3(dcmotW->min(0), muslinW->min(0), musfibW->min(0));
  double p_max = MAX3(dcmotW->max(0), muslinW->max(0), musfibW->max(0));

  double v_min = MIN3(dcmotW->min(1), muslinW->min(1), musfibW->min(1));
  double v_max = MAX3(dcmotW->max(1), muslinW->max(1), musfibW->max(1));

  double a_min = MIN3(dcmotW->min(2), muslinW->min(2), musfibW->min(2));
  double a_max = MAX3(dcmotW->max(2), muslinW->max(2), musfibW->max(2));

  double ac_min = dcmotA->min(0);
  double ac_max = dcmotA->max(0);

  cout << "Domains:" << endl;
  cout << "  Position:          " << p_min  << " " << p_max  << endl;
  cout << "  Velocity:          " << v_min  << " " << v_max  << endl;
  cout << "  Acceleration:      " << a_min  << " " << a_max  << endl;
  cout << "  Action (DC Motor): " << ac_min << " " << ac_max << endl;
  cout << "  Bins:              " << BINS   << endl;

  // normalising DCMot
  dcmotW->normaliseColumn(0, p_min,  p_max);
  dcmotW->normaliseColumn(1, v_min,  v_max);
  dcmotW->normaliseColumn(2, a_min,  a_max);

  dcmotA->normaliseColumn(0, ac_min, ac_max);

  // normalising MusLin
  muslinW->normaliseColumn(0, p_min,  p_max);
  muslinW->normaliseColumn(1, v_min,  v_max);
  muslinW->normaliseColumn(2, a_min,  a_max);

  // normalising MusFib
  musfibW->normaliseColumn(0, p_min,  p_max);
  musfibW->normaliseColumn(1, v_min,  v_max);
  musfibW->normaliseColumn(2, a_min,  a_max);

  double **wdomain = new double*[3];
  for(int i = 0; i < 3; i++)
  {
    wdomain[i] = new double[2];
    wdomain[i][0] = 0.0;
    wdomain[i][1] = 1.0;
  }
  int *wbins = new int[3];
  wbins[0]   = BINS;
  wbins[1]   = BINS;
  wbins[2]   = BINS;

  double **adomain = new double*[1];
  for(int i = 0; i < 1; i++)
  {
    adomain[i] = new double[2];
    adomain[i][0] = 0.0;
    adomain[i][1] = 1.0;
  }
  int *abins = new int[1];
  abins[0]   = BINS;

  dcmotW->setDomains(wdomain);
  dcmotW->setBinSizes(wbins);
  dcmotA->setDomains(adomain);
  dcmotA->setBinSizes(abins);

  muslinW->setDomains(wdomain);
  muslinW->setBinSizes(wbins);
  muslinA->setDomains(adomain);
  muslinA->setBinSizes(abins);

  musfibW->setDomains(wdomain);
  musfibW->setBinSizes(wbins);
  musfibA->setDomains(adomain);
  musfibA->setBinSizes(abins);

  cout << "discretising dcmot W" << endl;
  ULContainer *DdcmotW  = dcmotW->discretiseByColumn(); delete dcmotW;
  cout << "discretising dcmot A" << endl;
  ULContainer *DdcmotA  = dcmotA->discretiseByColumn(); delete dcmotA;
  ULContainer *dcW2 = DdcmotW->drop(1);
  ULContainer *dcW1 = DdcmotW->drop(-1);
  ULContainer *dcA1 = DdcmotA->drop(-1);

  cout << "Calculating DCMot:" << endl;
  double     dcmc_w     = sparse::MC_W(dcW2, dcW1, dcA1);
  cout << "  " << dcmc_w << endl;

  ULContainer *dcZ  = dcW2;
  ULContainer *dcXY = dcW1;
  (*dcXY) += *dcA1; // 2 columns. 1st column = A1, 2nd column = W1.

  vector<vector<int> > px;
  vector<vector<int> > py;

  vector<int> aa;
  aa.push_back(0); // W
  aa.push_back(1); // A
  px.push_back(aa);

  vector<int> bb;
  bb.push_back(0); // W'
  py.push_back(bb);

  vector<vector<int> > qx;
  vector<vector<int> > qy;

  vector<int> qa;
  qa.push_back(1); // A
  qx.push_back(qa);

  vector<int> qb;
  qb.push_back(0); // W'
  qy.push_back(qb);

  vector<Feature*> pfeatures;
  pfeatures.push_back(new Feature(0,0));

  vector<Feature*> qfeatures;
  qfeatures.push_back(new Feature(0,0));

  GIS* p = new GIS();
  p->setData(dcXY, dcZ);
  p->setFeatures(px,py,pfeatures);
  cout << "p init" << endl;
  p->init();

  GIS* q = new GIS();
  q->setData(dcXY, dcZ);
  q->setFeatures(qx,qy,qfeatures);
  cout << "q init" << endl;
  q->init();

  p->iterate();
  q->iterate();

  for(int i = 0; i < 500; i++)
  {
    if(p->error() > ERROR_THRESHOLD) p->iterate();
    if(q->error() > ERROR_THRESHOLD) q->iterate();
     if(i % 20 == 0)
    {
      cout << "p error (" << i << "): " << p->error() << endl;
      cout << "p: " << *p << endl;
      cout << "q error (" << i << "): " << q->error() << endl;
      cout << "q: " << *q << endl;
      // KL* kl = new KL(p, q);
      // cout << "kl: "  << kl->divergence2() << endl;
      // delete kl;
    }
    if(i % 20 == 0 && i > 0)
    {
      KL* kl = new KL(p, q);
      cout << "after " << i << " iterations: " << kl->divergence2() << endl;
      delete kl;
    }
  }

  // for(int i = 0; i < 50000; i++)
  // {
    // p->iterate();
    // if(i % 10 == 0) cout << "p error (" << i << "): " << p->error() << endl;
    // if(p->error() < ERROR_THRESHOLD) break;
  // }
  // cout << "p converged" << endl;

  // for(int i = 0; i < 50000; i++)
  // {
    // q->iterate();
    // if(i % 10 == 0) cout << "q error (" << i << "): " << q->error() << endl;
    // if(q->error() < ERROR_THRESHOLD) break;
  // }
  // cout << "q converged" << endl;

  KL* kl = new KL(p, q);

  cout << kl->divergence2() << endl;
}

