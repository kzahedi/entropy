#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE scgis_test
#include <boost/test/unit_test.hpp>

#include <entropy++/Container.h>
#include <entropy++/Csv.h>
#include <entropy++/iterativescaling/Model.h>
#include <entropy++/iterativescaling/SCGIS.h>
#include <entropy++/iterativescaling/KL.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/SparseMatrix.h>

#include <fstream>
#include <iostream>
#include <string>

#include <math.h>

# define EPSILON         0.01
# define ERROR_THRESHOLD 0.01

# define BINS 300
// # define BINS 10


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

BOOST_AUTO_TEST_SUITE(MorphologicalComputation)

BOOST_AUTO_TEST_CASE(MC_W)
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
  ULContainer *DdcmotW  = dcmotW->discretise(); delete dcmotW;
  cout << "discretising dcmot A" << endl;
  ULContainer *DdcmotA  = dcmotA->discretise(); delete dcmotA;
  ULContainer *dcW2 = DdcmotW->drop(1);
  ULContainer *dcW1 = DdcmotW->drop(-1);
  ULContainer *dcA1 = DdcmotA->drop(-1);

  cout << "Calculating DCMot:" << endl;
  double dcmc_w = sparse::MC_W(dcW2, dcW1, dcA1);
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

  SCGIS* p = new SCGIS();
  p->setData(dcXY, dcZ);
  p->setFeatures(px,py,pfeatures);
  cout << "p init" << endl;
  p->init();

  SCGIS* q = new SCGIS();
  q->setData(dcXY, dcZ);
  q->setFeatures(qx,qy,qfeatures);
  cout << "q init" << endl;
  q->init();

  p->iterate();
  q->iterate();

  for(int i = 0; i < 20000; i++)
  {
    if(p->error() > ERROR_THRESHOLD) p->iterate();
    if(q->error() > ERROR_THRESHOLD) q->iterate();
    // cout << "Iteration " << i << " error: " << p->error() << " : " << q->error() << endl;
    if(i % 100 == 0 && i > 0)
    {
      KL* kl = new KL(p, q);
      cout << "after " << i << " iterations: " << kl->divergence2() << endl;
      delete kl;
    }
  }

  // KL* kl = new KL(p, q);
  // cout << kl->divergence2() << endl;
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(Logic)
BOOST_AUTO_TEST_CASE(AND)
{
  ULContainer *xData = new ULContainer(4,2);
  *xData << 0 << 0;
  *xData << 0 << 1;
  *xData << 1 << 0;
  *xData << 1 << 1;

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

  SCGIS* independentModel = new SCGIS();
  independentModel->setData(xData, yData);
  independentModel->setFeatures(ia,ib,features);
  independentModel->init();

  for(int i = 0; i < 20000; i++)
  {
    independentModel->iterate();
    if(independentModel->error() < 0.000000001) break;
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

  for(int i = 0; i < 20000; i++)
  {
    dependentModel->iterate();
    if(dependentModel->error() < 0.000000001) break;
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
  cout << "AND (bits): " << kl->divergence2() << endl;
  cout << "AND (nats): " << kl->divergenceN() << endl;
}

BOOST_AUTO_TEST_CASE(ANDCMI)
{
  ULContainer *xData = new ULContainer(4,2);
  *xData << 0 << 0;
  *xData << 0 << 1;
  *xData << 1 << 0;
  *xData << 1 << 1;

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
  iaa.push_back(1);
  ia.push_back(iaa);

  vector<int> ibb;
  ibb.push_back(0);
  ib.push_back(ibb);

  vector<Feature*> features;
  features.push_back(new Feature(0,0));

  SCGIS* independentModel = new SCGIS();
  independentModel->setData(xData, yData);
  independentModel->setFeatures(ia,ib,features);
  independentModel->init();

  for(int i = 0; i < 5000; i++)
  {
    independentModel->iterate();
    if(independentModel->error() < 0.000000001) break;
  }

  independentModel->calculateProbabilities();

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

  for(int i = 0; i < 5000; i++)
  {
    dependentModel->iterate();
    if(dependentModel->error() < 0.000000001) break;
  }

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

  dependentModel->calculateProbabilities();
  BOOST_CHECK((fabs(dpycx(0,0) - dependentModel->p_y_c_x(0,0))   < EPSILON));
  BOOST_CHECK((fabs(dpycx(0,1) - dependentModel->p_y_c_x(0,1))   < EPSILON));
  BOOST_CHECK((fabs(dpycx(0,2) - dependentModel->p_y_c_x(0,2))   < EPSILON));
  BOOST_CHECK((fabs(dpycx(0,3) - dependentModel->p_y_c_x(0,3))   < EPSILON));
  BOOST_CHECK((fabs(dpycx(1,0) - dependentModel->p_y_c_x(1,0))   < EPSILON));
  BOOST_CHECK((fabs(dpycx(1,1) - dependentModel->p_y_c_x(1,1))   < EPSILON));
  BOOST_CHECK((fabs(dpycx(1,2) - dependentModel->p_y_c_x(1,2))   < EPSILON));
  BOOST_CHECK((fabs(dpycx(1,3) - dependentModel->p_y_c_x(1,3))   < EPSILON));
  BOOST_CHECK((fabs(dpx(0,0)   - dependentModel->p_x(0))         < EPSILON));
  BOOST_CHECK((fabs(dpx(0,1)   - dependentModel->p_x(1))         < EPSILON));
  BOOST_CHECK((fabs(dpx(0,2)   - dependentModel->p_x(2))         < EPSILON));
  BOOST_CHECK((fabs(dpx(0,3)   - dependentModel->p_x(3))         < EPSILON));

  BOOST_CHECK((fabs(ipycx(0,0) - independentModel->p_y_c_x(0,0)) < EPSILON));
  BOOST_CHECK((fabs(ipycx(0,1) - independentModel->p_y_c_x(0,1)) < EPSILON));
  BOOST_CHECK((fabs(ipycx(1,0) - independentModel->p_y_c_x(1,0)) < EPSILON));
  BOOST_CHECK((fabs(ipycx(1,1) - independentModel->p_y_c_x(1,1)) < EPSILON));
  BOOST_CHECK((fabs(ipx(0,0)   - independentModel->p_x(0))       < EPSILON));
  BOOST_CHECK((fabs(ipx(0,1)   - independentModel->p_x(1))       < EPSILON));

  ////////////////////////////////////////////////////////////////////////////////
  // Synergy
  ////////////////////////////////////////////////////////////////////////////////

  KL* kl = new KL(dependentModel, independentModel);
  cout << "AND CMI (bits): " << kl->divergence2() << endl;
  cout << "AND CMI (nats): " << kl->divergenceN() << endl;
}

BOOST_AUTO_TEST_CASE(OR)
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

  SCGIS* independentModel = new SCGIS();
  independentModel->setData(xData, yData);
  independentModel->setFeatures(ia,ib,features);
  independentModel->init();

  for(int i = 0; i < 20000; i++)
  {
    independentModel->iterate();
    if(independentModel->error() < ERROR_THRESHOLD) break;
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

  for(int i = 0; i < 20000; i++)
  {
    dependentModel->iterate();
    if(dependentModel->error() < ERROR_THRESHOLD) break;
  }

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

  dependentModel->calculateProbabilities();
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
  cout << "OR (bits): " << kl->divergence2() << endl;
  cout << "OR (nats): " << kl->divergenceN() << endl;
}

BOOST_AUTO_TEST_CASE(XOR)
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

  vector<int> ibb;
  ibb.push_back(0);
  ib.push_back(ibb);

  vector<Feature*> features;
  features.push_back(new Feature(0,0));
  features.push_back(new Feature(1,0));

  SCGIS* independentModel = new SCGIS();
  independentModel->setData(xData, yData);
  independentModel->setFeatures(ia,ib,features);
  independentModel->init();

  for(int i = 0; i < 10; i++)
  {
    independentModel->iterate();
    if(independentModel->error() < ERROR_THRESHOLD) break;
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

  for(int i = 0; i < 20000; i++)
  {
    dependentModel->iterate();
    independentModel->iterate();
  }


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

  dependentModel->calculateProbabilities();
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
  // BOOST_CHECK_CLOSE(EPSILON, kl->divergence(), EPSILON);
  cout << "XOR (bits): " << kl->divergence2() << endl;
  cout << "XOR (nats): " << kl->divergenceN() << endl;
}
BOOST_AUTO_TEST_SUITE_END()


