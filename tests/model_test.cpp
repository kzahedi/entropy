#include "model_test.h"

#include <entropy++/Container.h>
#include <entropy++/Csv.h>
#include <entropy++/iterativescaling/Model.h>
#include <entropy++/iterativescaling/GIS.h>

#include <fstream>
#include <iostream>
#include <string>

#include <math.h>

using namespace std;
using namespace entropy;

#ifndef PARENT
#define PARENT "/Users/zahedi/projects/entropy/experiments/hopping/data/"
#endif

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

using namespace std;
using namespace entropy::iterativescaling;


// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( modelTest );

void modelTest::testUnique()
{
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

  int bins = 300;
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

  GIS* dcmodel = new GIS();
  dcmodel->setData(dcXY, dcZ);
  dcmodel->setFeatures(a,b,features);
  dcmodel->createUniqueContainer();
  dcmodel->countObservedFeatures();

  CPPUNIT_ASSERT(dcmodel->X()->rows() > 0);
  CPPUNIT_ASSERT(dcmodel->X()->rows() == dcXY->rows());
  CPPUNIT_ASSERT(dcmodel->Y()->rows() > 0);
  CPPUNIT_ASSERT(dcmodel->Y()->rows() == dcZ->rows());

  CPPUNIT_ASSERT(dcmodel->uniqueX(0)->rows() < dcmodel->X()->rows());
  CPPUNIT_ASSERT(dcmodel->uniqueY(0)->rows() < dcmodel->Y()->rows());

  // for(int i = 0; i < a.size(); i++)
  // {
    // cout << "X unique: " << i << endl << *(dcmodel->uniqueX(i)) << endl;
  // }

  // for(int i = 0; i < b.size(); i++)
  // {
    // cout << "Y unique: " << i << endl << *(dcmodel->uniqueY(i)) << endl;
  // }

  cout << "Nr of features: " << dcmodel->nrOfFeatures() << endl;
    // for(int i = 0; i < dcmodel->nrOfFeatures(); i++)
    // {
      // Feature *f = (dcmodel->feature(i));
      // cout << "Feature " << i << ":" << endl << *f << endl;
    // }

  for(int i = 0; i < 500; i++)
  {
    dcmodel->iterate();
    if(dcmodel->error() < 0.00000001)
    {
      cout << "Error: " << dcmodel->error() << " with " << i << " iterations" << endl;
      break;
    }
  }

  // for(int i = 0; i < dcmodel->nrOfFeatures(); i++)
  // {
    // cout << (*dcmodel->feature(i)) << endl;
  // }

  ////////////////////////////////////////////////////////////////////////////////
  // MusFib MOTOR
  ////////////////////////////////////////////////////////////////////////////////
  ULContainer *mfZ  = mfW2;
  ULContainer *mfXY = mfW1;

  GIS* mfmodel = new GIS();
  mfmodel->setData(mfXY, mfZ);
  mfmodel->setFeatures(a,b,features);
  mfmodel->createUniqueContainer();
  mfmodel->countObservedFeatures();

  CPPUNIT_ASSERT(mfmodel->X()->rows() > 0);
  CPPUNIT_ASSERT(mfmodel->X()->rows() == mfXY->rows());
  CPPUNIT_ASSERT(mfmodel->Y()->rows() > 0);
  CPPUNIT_ASSERT(mfmodel->Y()->rows() == mfZ->rows());

  CPPUNIT_ASSERT(mfmodel->uniqueX(0)->rows() < mfmodel->X()->rows());
  CPPUNIT_ASSERT(mfmodel->uniqueY(0)->rows() < mfmodel->Y()->rows());

  for(int i = 0; i < 500; i++)
  {
    mfmodel->iterate();
    if(mfmodel->error() < 0.00000001)
    {
      cout << "Error: " << mfmodel->error() << " with " << i << " iterations" << endl;
      break;
    }
  }

  // cout << "Nr of features: " << mfmodel->nrOfFeatures() << endl;

  // for(int i = 0; i < mfmodel->nrOfFeatures(); i++)
  // {
    // cout << (*mfmodel->feature(i)) << endl;
  // }

  ////////////////////////////////////////////////////////////////////////////////
  // MusLin MOTOR
  ////////////////////////////////////////////////////////////////////////////////
  ULContainer *mlZ  = mlW2;
  ULContainer *mlXY = mlW1;

  GIS* mlmodel = new GIS();
  mlmodel->setData(mlXY, mlZ);
  mlmodel->setFeatures(a,b,features);
  mlmodel->createUniqueContainer();
  mlmodel->countObservedFeatures();

  CPPUNIT_ASSERT(mlmodel->X()->rows() > 0);
  CPPUNIT_ASSERT(mlmodel->X()->rows() == mlXY->rows());
  CPPUNIT_ASSERT(mlmodel->Y()->rows() > 0);
  CPPUNIT_ASSERT(mlmodel->Y()->rows() == mlZ->rows());

  CPPUNIT_ASSERT(mlmodel->uniqueX(0)->rows() < mlmodel->X()->rows());
  CPPUNIT_ASSERT(mlmodel->uniqueY(0)->rows() < mlmodel->Y()->rows());

  cout << "Nr of features: " << mlmodel->nrOfFeatures() << endl;

  for(int i = 0; i < 500; i++)
  {
    mlmodel->iterate();
    if(mlmodel->error() < 0.00000001)
    {
      cout << "Error: " << mlmodel->error() << " with " << i << " iterations" << endl;
      break;
    }
  }

}

void modelTest::testUnique2()
{
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

  GIS* dcmodel = new GIS();
  dcmodel->setData(dcXY, dcZ);
  dcmodel->setFeatures(a,b,features);
  dcmodel->createUniqueContainer();
  dcmodel->countObservedFeatures();

  CPPUNIT_ASSERT(dcmodel->X()->rows() > 0);
  CPPUNIT_ASSERT(dcmodel->X()->rows() == dcXY->rows());
  CPPUNIT_ASSERT(dcmodel->Y()->rows() > 0);
  CPPUNIT_ASSERT(dcmodel->Y()->rows() == dcZ->rows());

  CPPUNIT_ASSERT(dcmodel->uniqueX(0)->rows() < dcmodel->X()->rows());
  CPPUNIT_ASSERT(dcmodel->uniqueY(0)->rows() < dcmodel->Y()->rows());

  cout << "Nr of features: " << dcmodel->nrOfFeatures() << endl;

  for(int i = 0; i < 500; i++)
  {
    dcmodel->iterate();
    if(dcmodel->error() < 0.00000001)
    {
      cout << "Error: " << dcmodel->error() << " with " << i << " iterations" << endl;
      break;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////
  // MusFib MOTOR
  ////////////////////////////////////////////////////////////////////////////////
  ULContainer *mfZ  = mfW2;
  ULContainer *mfXY = mfW1;

  GIS* mfmodel = new GIS();
  mfmodel->setData(mfXY, mfZ);
  mfmodel->setFeatures(a,b,features);
  mfmodel->createUniqueContainer();
  mfmodel->countObservedFeatures();

  CPPUNIT_ASSERT(mfmodel->X()->rows() > 0);
  CPPUNIT_ASSERT(mfmodel->X()->rows() == mfXY->rows());
  CPPUNIT_ASSERT(mfmodel->Y()->rows() > 0);
  CPPUNIT_ASSERT(mfmodel->Y()->rows() == mfZ->rows());

  CPPUNIT_ASSERT(mfmodel->uniqueX(0)->rows() < mfmodel->X()->rows());
  CPPUNIT_ASSERT(mfmodel->uniqueY(0)->rows() < mfmodel->Y()->rows());

  cout << "Nr of features: " << mfmodel->nrOfFeatures() << endl;

  for(int i = 0; i < 500; i++)
  {
    mfmodel->iterate();
    if(mfmodel->error() < 0.00000001)
    {
      cout << "Error: " << mfmodel->error() << " with " << i << " iterations" << endl;
      break;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////
  // MusLin MOTOR
  ////////////////////////////////////////////////////////////////////////////////
  ULContainer *mlZ  = mlW2;
  ULContainer *mlXY = mlW1;

  GIS* mlmodel = new GIS();
  mlmodel->setData(mlXY, mlZ);
  mlmodel->setFeatures(a,b,features);
  mlmodel->createUniqueContainer();
  mlmodel->countObservedFeatures();
  mlmodel->generateExpected();

  CPPUNIT_ASSERT(mlmodel->X()->rows() > 0);
  CPPUNIT_ASSERT(mlmodel->X()->rows() == mlXY->rows());
  CPPUNIT_ASSERT(mlmodel->Y()->rows() > 0);
  CPPUNIT_ASSERT(mlmodel->Y()->rows() == mlZ->rows());

  CPPUNIT_ASSERT(mlmodel->uniqueX(0)->rows() < mlmodel->X()->rows());
  CPPUNIT_ASSERT(mlmodel->uniqueY(0)->rows() < mlmodel->Y()->rows());

  cout << "Nr of features: " << mlmodel->nrOfFeatures() << endl;

  for(int i = 0; i < 500; i++)
  {
    mlmodel->iterate();
    if(mlmodel->error() < 0.00000001)
    {
      cout << "Error: " << mlmodel->error() << " with " << i << " iterations" << endl;
      break;
    }
  }

}
