#include <entropy++/Csv.h>
#include <entropy++/sparse/MC_MI.h>
#include <entropy++/sparse/MC_W.h>
#include <entropy++/sparse/state/MC_W.h>
#include <entropy++/sparse/state/MC_MI.h>

#include <string>
#include <iostream>
#include <math.h>

#include <boost/filesystem.hpp>

using namespace std;
using namespace boost::filesystem;
using namespace entropy::sparse;
using namespace entropy::sparse::state;

# define MAX(a,b)    ((a>b)?a:b)
# define MAX3(a,b,c) MAX(a,MAX(b,c))

# define MIN(a,b)    ((a<b)?a:b)
# define MIN3(a,b,c) MIN(a,MIN(b,c))

# define POS 0
# define VEL 1
# define ACC 2
# define ACT 3
# define MSI 4 // muscle sensor input

// position      = SimData(:,2);
// velocity      = SimData(:,3);
// accelaration  = SimData(:,4);
// action        = SimData(:,10);
// muscle_sensor = SimData(:,5);

# define POS_MATLAB_INDEX                 (2  - 1)
# define VEL_MATLAB_INDEX                 (3  - 1)
# define ACC_MATLAB_INDEX                 (4  - 1)
# define ACT_MATLAB_INDEX                 (10 - 1)
# define MUSCLE_SENSOR_INPUT_MATLAB_INDEX (5  - 1)

int main(int argc, char **argv)
{
  if(argc != 2)
  {
    cout << "usage: " << argv[0] << " <directory>" << endl;
    exit(-1);
  }

  string d = string(argv[1]);

  Csv *csv = new Csv();

  Container *dcmot  = csv->read(d + "/dcmot.csv",  4,
                                POS_MATLAB_INDEX,
                                VEL_MATLAB_INDEX,
                                ACC_MATLAB_INDEX,
                                ACT_MATLAB_INDEX);

  Container *muslin = csv->read(d + "/muslin.csv", 5,
                                POS_MATLAB_INDEX,
                                VEL_MATLAB_INDEX,
                                ACC_MATLAB_INDEX,
                                ACT_MATLAB_INDEX,
                                MUSCLE_SENSOR_INPUT_MATLAB_INDEX);

  Container *musfib = csv->read(d + "/musfib.csv", 5,
                                POS_MATLAB_INDEX,
                                VEL_MATLAB_INDEX,
                                ACC_MATLAB_INDEX,
                                ACT_MATLAB_INDEX,
                                MUSCLE_SENSOR_INPUT_MATLAB_INDEX);

  Container *dcmotW  = dcmot->columns(3, 0, 1, 2);
  Container *dcmotS  = dcmot->columns(2, 0, 1);
  Container *dcmotA  = dcmot->columns(1, 3);
  delete dcmot;

  Container *muslinW = muslin->columns(3, 0, 1, 2);
  Container *muslinS = muslin->columns(1, 4);
  Container *muslinA = muslin->columns(1, 3);
  delete muslin;

  Container *musfibW = musfib->columns(3, 0, 1, 2);
  Container *musfibS = musfib->columns(1, 4);
  Container *musfibA = musfib->columns(1, 3);
  delete musfib;

  double p_min = MIN3(dcmotW->min(0), muslinW->min(0), musfibW->min(0));
  double p_max = MAX3(dcmotW->max(0), muslinW->max(0), musfibW->max(0));

  double v_min = MIN3(dcmotW->min(1), muslinW->min(1), musfibW->min(1));
  double v_max = MAX3(dcmotW->max(1), muslinW->max(1), musfibW->max(1));

  double a_min = MIN3(dcmotW->min(2), muslinW->min(2), musfibW->min(2));
  double a_max = MAX3(dcmotW->max(2), muslinW->max(2), musfibW->max(2));

  double mi_min = MIN(muslinS->min(0), musfibS->min(0));
  double mi_max = MAX(muslinS->max(0), musfibS->max(0));

  double ac_min = dcmotA->min(2);
  double ac_max = dcmotA->max(2);

  cout << "Domains:" << endl;
  cout << "  Position:          " << p_min  << " " << p_max  << endl;
  cout << "  Velocity:          " << v_min  << " " << v_max  << endl;
  cout << "  Acceleration:      " << a_min  << " " << a_max  << endl;
  cout << "  Action (DC Motor): " << ac_min << " " << ac_max << endl;

  // normalising DCMot
  dcmotW->normaliseColumn(0, p_min,  p_max);
  dcmotW->normaliseColumn(1, v_min,  v_max);
  dcmotW->normaliseColumn(2, a_min,  a_max);

  dcmotS->normaliseColumn(0, p_min,  p_max);
  dcmotS->normaliseColumn(1, v_min,  v_max);

  dcmotA->normaliseColumn(0, ac_min, ac_max);

  // normalising MusLin
  muslinW->normaliseColumn(0, p_min,  p_max);
  muslinW->normaliseColumn(1, v_min,  v_max);
  muslinW->normaliseColumn(2, a_min,  a_max);

  muslinS->normaliseColumn(0, mi_min, mi_max);

  // normalising MusFib
  musfibW->normaliseColumn(0, p_min,  p_max);
  musfibW->normaliseColumn(1, v_min,  v_max);
  musfibW->normaliseColumn(2, a_min,  a_max);

  musfibS->normaliseColumn(0, mi_min, mi_max);

  double **wdomain = new double*[3];
  for(int i = 0; i < 3; i++)
  {
    wdomain[i] = new double[2];
    wdomain[i][0] = 0.0;
    wdomain[i][1] = 1.0;
  }
  int bins   = 300;
  int *wbins = new int[3];
  wbins[0]   = bins;
  wbins[1]   = bins;
  wbins[2]   = bins;

  double **sdomain = new double*[2];
  for(int i = 0; i < 2; i++)
  {
    sdomain[i] = new double[2];
    sdomain[i][0] = 0.0;
    sdomain[i][1] = 1.0;
  }
  int *sbins = new int[2];
  sbins[0]   = bins;
  sbins[1]   = bins;

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
  dcmotS->setDomains(sdomain);
  dcmotS->setBinSizes(sbins);
  dcmotA->setDomains(adomain);
  dcmotA->setBinSizes(abins);

  muslinW->setDomains(wdomain);
  muslinW->setBinSizes(wbins);
  muslinS->setDomains(sdomain);
  muslinS->setBinSizes(sbins);
  muslinA->setDomains(adomain);
  muslinA->setBinSizes(abins);

  musfibW->setDomains(wdomain);
  musfibW->setBinSizes(wbins);
  musfibS->setDomains(sdomain);
  musfibS->setBinSizes(sbins);
  musfibA->setDomains(adomain);
  musfibA->setBinSizes(abins);

  cout << "discretising dcmot W" << endl;
  Container *DdcmotW  = dcmotW->discretise(); delete dcmotW;
  cout << "discretising dcmot S" << endl;
  Container *DdcmotS  = dcmotS->discretise(); delete dcmotS;
  cout << "discretising dcmot A" << endl;
  Container *DdcmotA  = dcmotA->discretise(); delete dcmotA;
  Container *dcW2 = DdcmotW->drop(1);
  Container *dcW1 = DdcmotW->drop(-1);
  Container *dcA1 = DdcmotA->drop(-1);
  Container *dcS1 = DdcmotS->drop(-1);

  cout << "discretising muslin W" << endl;
  Container *DmuslinW  = muslinW->discretise(); delete muslinW;
  cout << "discretising muslin S" << endl;
  Container *DmuslinS  = muslinS->discretise(); delete muslinS;
  cout << "discretising muslin A" << endl;
  Container *DmuslinA  = muslinA->discretise(); delete muslinA;
  Container *mlW2 = DmuslinW->drop(1);
  Container *mlW1 = DmuslinW->drop(-1);
  Container *mlA1 = DmuslinA->drop(-1);
  Container *mlS1 = DmuslinS->drop(-1);

  cout << "discretising musfib W" << endl;
  Container *DmusfibW  = musfibW->discretise(); delete musfibW;
  cout << "discretising musfib S" << endl;
  Container *DmusfibS  = musfibS->discretise(); delete musfibS;
  cout << "discretising musfib A" << endl;
  Container *DmusfibA  = musfibA->discretise(); delete musfibA;
  Container *mfW2 = DmusfibW->drop(1);
  Container *mfW1 = DmusfibW->drop(-1);
  Container *mfA1 = DmusfibA->drop(-1);
  Container *mfS1 = DmusfibS->drop(-1);

  cout << "Calculating DCMot:" << endl;
  double     dcmc_w   = entropy::sparse::MC_W(dcW2, dcW1, dcA1);
  double     dcmc_mi  = entropy::sparse::MC_MI(dcW2, dcW1, dcS1, dcA1);
  Container* dcmc_mid = entropy::sparse::state::MC_MI(dcW2, dcW1, dcS1, dcA1);
  Container* dcmc_wd  = entropy::sparse::state::MC_W(dcW2, dcW1, dcA1);
  double     s        = dcmc_mid->columnSum(0) / ((double)dcmc_mid->rows());
  double     ss       = dcmc_wd->columnSum(0)  / ((double)dcmc_mid->rows());
  cout << "  DCMot MC_W:  " << dcmc_w  << endl;
  cout << "  DCMot MC_MI: " << dcmc_mi << endl;
  cout << "  checking values: " << s  << " vs. " << dcmc_mi << " = " << fabs(s  - dcmc_mi) << endl;
  cout << "  checking values: " << ss << " vs. " << dcmc_w  << " = " << fabs(ss - dcmc_w)  << endl;
  cout << "writing DCMot state-dependent MC_W" << endl;
  csv->write("state-dependent-dc_mc_w.csv", dcmc_wd);
  cout << "writing DCMot state-dependent MC_MI" << endl;
  csv->write("state-dependent-dc_mc_mi.csv", dcmc_mid);

  cout << "Calculating MusLin:" << endl;
  double     mlmc_w   = entropy::sparse::MC_W(mlW2, mlW1, mlA1);
  double     mlmc_mi  = entropy::sparse::MC_MI(mlW2, mlW1, mlS1, mlA1);
  Container* mlmc_mid = entropy::sparse::state::MC_MI(mlW2, mlW1, mlS1, mlA1);
  Container* mlmc_wd  = entropy::sparse::state::MC_W(mlW2, mlW1, mlA1);
  double     t        = mlmc_mid->columnSum(0) / ((double)mlmc_mid->rows());
  double     tt       = mlmc_wd->columnSum(0) / ((double)mlmc_wd->rows());
  cout << "MusLin MC_W:  " << mlmc_w  << endl;
  cout << "MusLin MC_MI: " << mlmc_mi << endl;
  cout << "  checking values: " << t  << " vs. " << mlmc_mi << " = " << fabs(t  - mlmc_mi) << endl;
  cout << "  checking values: " << tt << " vs. " << mlmc_w  << " = " << fabs(tt - mlmc_w)  << endl;
  cout << "writing MusLin state-dependent MC_W" << endl;
  csv->write("state-dependent-ml_mc_w.csv", mlmc_wd);
  cout << "writing MusLin state-dependent MC_MI" << endl;
  csv->write("state-dependent-ml_mc_mi.csv", mlmc_mid);

  cout << "Calculating MusFib:" << endl;
  double     mfmc_w   = entropy::sparse::MC_W(mfW2, mfW1, mfA1);
  double     mfmc_mi  = entropy::sparse::MC_MI(mfW2, mfW1, mfS1, mfA1);
  Container* mfmc_mid = entropy::sparse::state::MC_MI(mfW2, mfW1, mfS1, mfA1);
  Container* mfmc_wd  = entropy::sparse::state::MC_W(mfW2, mfW1, mfA1);
  double     u        = mfmc_mid->columnSum(0) / ((double)mfmc_mid->rows());
  double     uu       = mfmc_wd->columnSum(0)  / ((double)mfmc_mid->rows());
  cout << "MusFib MC_W:  " << mfmc_w  << endl;
  cout << "MusFib MC_MI: " << mfmc_mi << endl;
  cout << "  checking values: " << u  << " vs. " << mfmc_mi << " = " << fabs(u  - mfmc_mi) << endl;
  cout << "  checking values: " << uu << " vs. " << mfmc_w  << " = " << fabs(uu - mfmc_w)  << endl;
  cout << "writing MusFib state-dependent MC_W" << endl;
  csv->write("state-dependent-mf_mc_w.csv", mfmc_wd);
  cout << "writing MusFib state-dependent MC_MI" << endl;
  csv->write("state-dependent-mf_mc_mi.csv", mfmc_mid);

}
