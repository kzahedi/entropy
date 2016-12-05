#include "it_cmi_test.h"

#include <entropy++/iterativescaling/gis/IterativeScaling.h>
#include <entropy++/iterativescaling/scgis/IterativeScaling.h>
// #include "../experiments/it/Test.h"
// #include "../experiments/it/TestMI.h"
#include <entropy++/Container.h>
#include <entropy++/Csv.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/sparse/state/MC_W.h>

#include <string>
#include <iostream>
#include <math.h>

#include <boost/filesystem.hpp>

#ifndef PARENT
#define PARENT "/Users/zahedi/projects/entropy/experiments/hopping/data/"
#endif

using namespace std;
using namespace boost::filesystem;
using namespace entropy::sparse;
using namespace entropy::sparse::state;
using namespace entropy::iterativescaling;

# define MAX3(a,b,c) MAX(a,MAX(b,c))
# define MIN3(a,b,c) MIN(a,MIN(b,c))

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

#include <iostream>
#include <string>

using namespace std;
using namespace entropy;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( itCmiTest );

void itCmiTest::testITvsCMI()
{
  cout << PARENT << "/dcmot.csv" << endl;

  Csv *csv = new Csv();
  DContainer *dcmot  = csv->read(string(PARENT) + string("/dcmot.csv"),  4,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX);

  DContainer *muslin = csv->read(string(PARENT) + string("/muslin.csv"), 5,
                                 POS_MATLAB_INDEX,
                                 VEL_MATLAB_INDEX,
                                 ACC_MATLAB_INDEX,
                                 ACT_MATLAB_INDEX,
                                 MUSCLE_SENSOR_INPUT_MATLAB_INDEX);

  DContainer *musfib = csv->read(string(PARENT) + string("/musfib.csv"), 5,
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

  int bins = 5;
  cout << "Domains:" << endl;
  cout << "  Position:          " << p_min  << " " << p_max  << endl;
  cout << "  Velocity:          " << v_min  << " " << v_max  << endl;
  cout << "  Acceleration:      " << a_min  << " " << a_max  << endl;
  cout << "  Action (DC Motor): " << ac_min << " " << ac_max << endl;
  cout << "  Bins:              " << bins   << endl;

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
  wbins[0]   = bins;
  wbins[1]   = bins;
  wbins[2]   = bins;

  double **adomain = new double*[1];
  for(int i = 0; i < 1; i++)
  {
    adomain[i] = new double[2];
    adomain[i][0] = 0.0;
    adomain[i][1] = 1.0;
  }
  int *abins = new int[1];
  abins[0]   = bins;

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

  cout << "discretising muslin W" << endl;
  ULContainer *DmuslinW  = muslinW->discretise(); delete muslinW;
  cout << "discretising muslin A" << endl;
  ULContainer *DmuslinA  = muslinA->discretise(); delete muslinA;
  ULContainer *mlW2 = DmuslinW->drop(1);
  ULContainer *mlW1 = DmuslinW->drop(-1);
  ULContainer *mlA1 = DmuslinA->drop(-1);

  cout << "discretising musfib W" << endl;
  ULContainer *DmusfibW  = musfibW->discretise(); delete musfibW;
  cout << "discretising musfib A" << endl;
  ULContainer *DmusfibA  = musfibA->discretise(); delete musfibA;
  ULContainer *mfW2 = DmusfibW->drop(1);
  ULContainer *mfW1 = DmusfibW->drop(-1);
  ULContainer *mfA1 = DmusfibA->drop(-1);

  cout << "Calculating DCMot:" << endl;
  double dcmc_w = entropy::sparse::MC_W(dcW2, dcW1, dcA1);
  cout << "  DCMot MC_W:  " << dcmc_w  << endl;

  cout << "Calculating MusLin:" << endl;
  double     mlmc_w   = entropy::sparse::MC_W(mlW2, mlW1, mlA1);
  cout << "MusLin MC_W:  " << mlmc_w  << endl;

  cout << "Calculating MusFib:" << endl;
  double     mfmc_w   = entropy::sparse::MC_W(mfW2, mfW1, mfA1);
  cout << "MusFib MC_W:  " << mfmc_w  << endl;

  // Tests start here

  // I(Y;X|Z) = I(W';W|A)
  ULContainer *dcZ  = dcW2;
  ULContainer *dcXY = dcW1;
  (*dcXY) += *dcA1; // 2 columns. 1st column = A1, 2nd column = W1. 
  // cout << (*dcZ) << endl;
  // cout << (*dcXY) << endl;
  CPPUNIT_ASSERT(dcW1->rows() > 0);
  for(int r = 0; r < dcW1->rows(); r++)
  {
    CPPUNIT_ASSERT_EQUAL((*dcXY)(r,0), (*dcW1)(r,0));
    CPPUNIT_ASSERT_EQUAL((*dcXY)(r,1), (*dcA1)(r,0));
  }

  vector<vector<int> > systAX(1,vector<int>(0));
  systAX[0].push_back(0);
  systAX[0].push_back(1); // x, y -> z

  vector<vector<int> > systAY(1,vector<int>(0));
  systAY[0].push_back(0); // z

  vector<vector<int> > systBX(1,vector<int>(0));
  systBX[0].push_back(1); // y -> z
  vector<vector<int> > systBY(1,vector<int>(0));
  systBY[0].push_back(0); // z

  ULContainer *xAlphabet = dcXY->globalUnique();
  ULContainer *yAlphabet = dcZ->globalUnique();

  // cout << *xAlphabet << endl;
  // cout << *yAlphabet << endl;

  IsParameter param;
  param.lambdavalue    = 1.0;
  param.lambdadeltaval = 1.0;
  param.sigma          = 0.01;   //TODO find best sigma value
  param.maxit          = 100;
  param.konv           = 0.0001;
  param.time           = false;
  param.test           = true;
  param.seconds        = 600;
  // param.verbose        = true;

  gis::IterativeScaling* p1 = new gis::IterativeScaling(dcXY, dcZ, xAlphabet, yAlphabet, systAX, systAY, param);
  gis::IterativeScaling* p2 = new gis::IterativeScaling(dcXY, dcZ, xAlphabet, yAlphabet, systBX, systBY, param);

  int I      = pow(xAlphabet->rows(), dcXY->columns());
  int J      = pow(yAlphabet->rows(), dcZ->columns());
  double p   = 0.0;
  double q   = 0.0;
  double val = 0.0;
  cout << "I: " << I << " J: " << J << endl;
  for(int i = 0; i < I; i++)
  {
    for( int j = 0; j < J; j++)
    {
      p = p1->p_x_c_y(i,j);
      q = p2->p_x_c_y(i,j);
      cout << q << endl;
      if(q > 0.000001)
      {
        double v = p * p1->propm(i) * log2(p / q);
        val += v;
        cout << v << " -> " << val << endl;
      }
    }
  }
  cout << val << " " << dcmc_w << endl;
}
