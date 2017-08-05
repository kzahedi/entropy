#include <entropy++/Random.h>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include <boost/filesystem.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/progress.hpp>
#include <boost/tokenizer.hpp>

#include <entropy++/Container.h>
#include <entropy++/Csv.h>
#include <entropy++/distribution/MC_A.h>
#include <entropy++/distribution/MC_MI.h>
#include <entropy++/distribution/MC_W.h>
#include <entropy++/iterativescaling/SCGIS.h>
#include <entropy++/iterativescaling/GIS.h>
#include <entropy++/iterativescaling/KL.h>
#include <entropy++/iterativescaling/Original.h>

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

#include <math.h>

#define LOOP(a,b) for(int a = 0; a < (int)b.size();a++)
#define BIN(a)    ((a == 0)?-1:1)
#define MIN(a,b)  ((a<b)?a:b)
#define MAX(a,b)  ((a>b)?a:b)
// #define MAP(a)    ((a<0)?0:((a>1)?1:a))
#define MAP(a)    a
#define NONE      -1
#define __MC_W        1001
#define __MC_A        1002
#define __MC_MI       1003
#define __MC_SY_SCGIS 1004
#define __MC_SY_GIS   1005
#define __MC_SYO      1006
#define __MC_SYO_NID  1007
#define __GENERATE    1008

using namespace std;
using namespace boost;
using namespace entropy;
using namespace entropy::distribution;
using namespace entropy::iterativescaling;

int    bins        = -1;
bool   logLambda   = false;
bool   logMatrices = false;

DEFINE_string(mu,   "0.0",      "mu values   (s -> a),    0.0 = random, 1.0 = deterministic");
DEFINE_string(phi,  "0.0",      "phi values  (w -> w'),   0.0 = random, 1.0 = deterministic");
DEFINE_string(psi,  "0.0",      "psi values  (a -> w'),   0.0 = random, 1.0 = deterministic");
DEFINE_string(chi,  "0.0",      "chi values  (a,w -> w'), 0.0 = random, 1.0 = deterministic");
DEFINE_string(zeta, "0.0",      "zeta values (w -> s),    0.0 = random, 1.0 = deterministic");
DEFINE_string(tau,  "0.0",      "tau values  (p(w)),      0.0 = random, 1.0 = deterministic");
DEFINE_bool(L,      false,      "Log the lambda values");
DEFINE_bool(M,      false,      "Log matrices");

DEFINE_string(mc,   "MC_W",     "quantification");
DEFINE_string(o,    "mc_w.csv", "output file");
DEFINE_int64(b,     2,          "bins per random variable");
DEFINE_int64(syci,  1000,       "SYO convergence iterations");
DEFINE_int64(sysi,  100,        "SY data drive sample iterations");


class Indicator
{
  public:
    Indicator(int a, int b)
    {
      _a      = a;
      _b      = b;
      _c      = -1;
      _lambda = 1.0;
    };

    Indicator(int a, int b, int c)
    {
      _a      = a;
      _b      = b;
      _c      = c;
      _lambda = 1.0;
    };

    int a() {return _a;};
    int b() {return _b;};
    int c() {return _c;};

    void setLambda(double v)
    {
      _lambda = v;
    }

    double lambda()
    {
      return _lambda;
    }

    double lambda(int a, int b)
    {
      if(a == _a && b == _b) return _lambda;
      return -_lambda;
    }

    double lambda(int a, int b, int c)
    {
      if(a == _a && b == _b && c == _c) return _lambda;
      return -_lambda;
    }

    friend std::ostream& operator<<(std::ostream& str, const Indicator& i)
    {
      str << "Indicator: ";
      str << i._a << ", " << i._b;
      if(i._c > -1) str << ", " << i._c;
      str << ": " << i._lambda;
      return str;
    };

  private:
    int    _a;
    int    _b;
    int    _c;
    double _lambda;
};

double get_value(double a, // double amax, double amin,
                 double b, // double bmin, double bmax,
                 double factor)
{
  double avalue = 2.0 * a/((float)(bins - 1.0)) - 1.0;
  double bvalue = 2.0 * b/((float)(bins - 1.0)) - 1.0;
  return factor * avalue * bvalue;
}

double get_value2(double a, // double amin, double amax,
                  double b, // double bmin, double bmax,
                  double c, // double cmin, double cmax,
                  double factor)
{
  double avalue = 2.0 * a/((float)(bins-1.0)) - 1.0;
  double bvalue = 2.0 * b/((float)(bins-1.0)) - 1.0;
  double cvalue = 2.0 * c/((float)(bins-1.0)) - 1.0;
  return factor * avalue * bvalue * cvalue;
}

vector<float> float_tokenizer(string s)
{
  vector<float> values;

  if(s.find(",") == std::string::npos)
  {
    boost::char_separator<char> sep(":");
    tokenizer<boost::char_separator<char> > tk(s,sep);
    vector<float> vec;
    for(tokenizer<boost::char_separator<char> >::iterator i(tk.begin()); i!=tk.end();++i)
    {
      vec.push_back(atof((*i).c_str()));
    }

    double start = vec[0];
    double delta = 0.01;
    double end   = vec[0];
    if(vec.size() > 1)
    {
      start = vec[0];
      delta = vec[1];
      end   = vec[2];
    }

    double value = start;

    while(value <= end)
    {
      values.push_back(value);
      value += delta;
    }
    if(values[values.size()-1] < end) values.push_back(end);
  }
  else
  {
    boost::char_separator<char> sep(",");
    tokenizer<boost::char_separator<char> > tk(s,sep);
    vector<float> vec;
    for(tokenizer<boost::char_separator<char> >::iterator i(tk.begin()); i!=tk.end();++i)
    {
      values.push_back(atof((*i).c_str()));
    }
  }
  return values;
}

void clear_pw1(double *a)
{
  if(a != NULL)
  {
    delete a;
  }
  a = NULL;
}

// double* create_pw1()
// {
  // double *a = new double[bins];
  // return a;
// }

void clear_pa1_c_s1(double **a)
{
  if(a != NULL)
  {
    for(int i = 0; i < bins; i++)
    {
      delete a[i];
    }
    delete a;
  }
  a = NULL;
}

// double** create_pa1_c_s1()
void create_pa1_c_s1(double** a)
{
  // double** a = new double*[bins];
  for(int i = 0; i < bins; i++)
  {
    a[i] = new double[bins];
    for(int j = 0; j < bins; j++)
    {
      a[i][j] = 0.0;
    }
  }
  // return a;
}

void clear_ps1_c_w1(double** a)
{
  if(a != NULL)
  {
    for(int i = 0; i < bins; i++)
    {
      delete a[i];
    }
    delete a;
  }
  a = NULL;
}

// double** create_ps1_c_w1()
void create_ps1_c_w1(double **a)
{
  // double **a = new double*[bins];
  for(int i = 0; i < bins; i++)
  {
    a[i] = new double[bins];
    for(int j = 0; j < bins; j++)
    {
      a[i][j] = 0.0;
    }
  }
  // return a;
}

void clear_pw2_c_w1a1(double*** a)
{
  if(a != NULL)
  {
    for(int i = 0; i < bins; i++)
    {
      for(int j = 0; j < bins; j++)
      {
        delete a[i][j];
      }
      delete a[i];
    }
    delete a;
  }
  a = NULL;
}

// double*** create_pw2_c_w1a1()
void create_pw2_c_w1a1(double*** a)
{
  // double*** a = new double**[bins];
  for(int i = 0; i < bins; i++)
  {
    a[i] = new double*[bins];
    for(int j = 0; j < bins; j++)
    {
      a[i][j] = new double[bins];
      for(int k = 0; k < bins; k++)
      {
        a[i][j][k] = 0.0;
      }
    }
  }
  // return a;
}

void clear_pw2w1s1a1(double**** a)
{
  if(a != NULL)
  {
    for(int i = 0; i < bins; i++)
    {
      for(int j = 0; j < bins; j++)
      {
        for(int k = 0; k < bins; k++)
        {
          delete a[i][j][k];
        }
        delete a[i][j];
      }
      delete a[i];
    }
    delete a;
  }
  a = NULL;
}

// double**** create_pw2w1s1a1(double**** a)
void create_pw2w1s1a1(double**** a)
{
  // double**** a = new double***[bins];
  // a = new double***[bins];
  for(int i = 0; i < bins; i++)
  {
    a[i] = new double**[bins];
    for(int j = 0; j < bins; j++)
    {
      a[i][j] = new double*[bins];
      for(int k = 0; k < bins; k++)
      {
        a[i][j][k] = new double[bins];
        for(int l = 0; l < bins; l++)
        {
          a[i][j][k][l] = 0.0;
        }
      }
    }
  }
  // return a;
}

// double*** create_pw2w1a1()
void create_pw2w1a1(double*** a)
{
  // double*** a = new double**[bins];
  for(int i = 0; i < bins; i++)
  {
    a[i] = new double*[bins];
    for(int j = 0; j < bins; j++)
    {
      a[i][j] = new double[bins];
      for(int k = 0; k < bins; k++)
      {
        a[i][j][k] = 0.0;
      }
    }
  }
  // return a;
}

void clear_pw2w1a1(double*** a)
{
  if(a != NULL)
  {
    for(int i = 0; i < bins; i++)
    {
      for(int j = 0; j < bins; j++)
      {
        delete a[i][j];
      }
      delete a[i];
    }
    delete a;
  }
  a = NULL;
}

// double*** create_pw2a1w1()
void create_pw2a1w1(double*** a)
{
  // double*** a = new double**[bins];
  for(int i = 0; i < bins; i++)
  {
    a[i] = new double*[bins];
    for(int j = 0; j < bins; j++)
    {
      a[i][j] = new double[bins];
      for(int k = 0; k < bins; k++)
      {
        a[i][j][k] = 0.0;
      }
    }
  }
  // return a;
}

void clear_pw2a1w1(double*** a)
{
  if(a != NULL)
  {
    for(int i = 0; i < bins; i++)
    {
      for(int j = 0; j < bins; j++)
      {
        delete a[i][j];
      }
      delete a[i];
    }
    delete a;
  }
  a = NULL;
}


double pw2_c_w1_a1(int w2, int w1, int a1,
                   vector<Indicator*>& i_w2w1a1,
                   vector<Indicator*>& i_w2w1,
                   vector<Indicator*>& i_w2a1)
{
  double z = 0.0;
  double n = 0.0;

  for(vector<Indicator*>::iterator i = i_w2w1a1.begin(); i != i_w2w1a1.end(); i++) z += (*i)->lambda(w2, w1, a1);
  for(vector<Indicator*>::iterator i = i_w2w1.begin();   i != i_w2w1.end();   i++) z += (*i)->lambda(w2, w1);
  for(vector<Indicator*>::iterator i = i_w2a1.begin();   i != i_w2a1.end();   i++) z += (*i)->lambda(w2, a1);

  for(int w22 = 0; w22 < bins; w22++)
  {
    double nn = 0.0;
    for(vector<Indicator*>::iterator i = i_w2w1a1.begin(); i != i_w2w1a1.end(); i++) nn += (*i)->lambda(w22, w1, a1);
    for(vector<Indicator*>::iterator i = i_w2w1.begin();   i != i_w2w1.end();   i++) nn += (*i)->lambda(w22, w1);
    for(vector<Indicator*>::iterator i = i_w2a1.begin();   i != i_w2a1.end();   i++) nn += (*i)->lambda(w22, a1);
    n += exp(nn);
  }

  return exp(z)/n;
}

double pa1_c_s1(int a1, int s1, vector<Indicator*>& i_a1s1)
{
  double z = 0.0;
  double n = 0.0;

  for(vector<Indicator*>::iterator i = i_a1s1.begin(); i != i_a1s1.end(); i++)
    z += (*i)->lambda(a1, s1);

  for(int a11 = 0; a11 < bins; a11++)
  {
    double nn = 0.0;
    for(vector<Indicator*>::iterator i = i_a1s1.begin(); i != i_a1s1.end(); i++)
      nn += (*i)->lambda(a11, s1);
    n += exp(nn);
  }

  return exp(z)/n;
}

double ps1_c_w1(int s1, int w1, vector<Indicator*>& i_s1w1)
{
  double z = 0.0;
  double n = 0.0;

  for(vector<Indicator*>::iterator i = i_s1w1.begin(); i != i_s1w1.end(); i++)
    z += (*i)->lambda(s1, w1);

  for(int s11 = 0; s11 < bins; s11++)
  {
    double nn = 0.0;
    for(vector<Indicator*>::iterator i = i_s1w1.begin(); i != i_s1w1.end(); i++)
      nn += (*i)->lambda(s11, w1);
    n += exp(nn);
  }
  return exp(z)/n;
}

double pw1(int w1)
{
  return 1.0/(float)(bins);
}

string toString(double**** pw2w1s1a1)
{
  stringstream sst;
  sst << "# w2, w1, s1, a1, value" << endl;
  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int s1 = 0; s1 < bins; s1++)
        for(int a1 = 0; a1 < bins; a1++)
          sst << BIN(w2) << "," << BIN(w1) << "," << BIN(s1) << "," << BIN(a1) << ","
            << pw2w1s1a1[w2][w1][s1][a1] << endl;
  return sst.str();
};

void generate_probability_distrbution(double mu,   double phi,
                                      double psi,  double chi,
                                      double zeta, double tau,
                                      double**** pw2w1s1a1,
                                      double***  m_pw2_c_w1_a1,
                                      double**   m_pa1_c_s1,
                                      double**   m_ps1_c_w1,
                                      double*    m_pw1)
{
  vector<Indicator*> i_w2w1a1;
  vector<Indicator*> i_w2w1;
  vector<Indicator*> i_w2a1;
  vector<Indicator*> i_s1w1;
  vector<Indicator*> i_a1s1;

  for(int w2 = 0; w2 < bins; w2++)
  {
    for(int w1 = 0; w1 < bins; w1++)
    {
      for(int a1 = 0; a1 < bins; a1++)
      {
        Indicator *i = new Indicator(w2, w1, a1);
        i->setLambda(get_value2(w2,w1,a1,chi));
        i_w2w1a1.push_back(i);
      }
    }
  }

  for(int w2 = 0; w2 < bins; w2++)
  {
    for(int w1 = 0; w1 < bins; w1++)
    {
      Indicator *i = new Indicator(w2, w1);
      i->setLambda(get_value(w2,w1,phi));
      i_w2w1.push_back(i);
    }
  }

  for(int w2 = 0; w2 < bins; w2++)
  {
    for(int a1 = 0; a1 < bins; a1++)
    {
      Indicator *i = new Indicator(w2, a1);
      i->setLambda(get_value(w2,a1,psi));
      i_w2a1.push_back(i);
    }
  }

  for(int a1 = 0; a1 < bins; a1++)
  {
    for(int s1 = 0; s1 < bins; s1++)
    {
      Indicator *i = new Indicator(a1, s1);
      i->setLambda(get_value(a1,s1,mu));
      i_a1s1.push_back(i);
    }
  }

  for(int s1 = 0; s1 < bins; s1++)
  {
    for(int w1 = 0; w1 < bins; w1++)
    {
      Indicator *i = new Indicator(s1, w1);
      i->setLambda(get_value(s1,w1,zeta));
      i_s1w1.push_back(i);
    }
  }

  // cout << "Indicators:" << endl;;
  if(logLambda)
  {
    stringstream sst;
    sst << "logs/phi_" << phi << "_psi_" << psi << "_chi_" << chi << "_mu_" << mu << "_zeta_" << zeta << "_tau_" << tau;
    string subdir = sst.str();
    create_directory(boost::filesystem::path(subdir));

    sst.str("");
    sst << subdir << "/" << "w2w1a1.csv";
    ofstream output(sst.str().c_str());
    for(vector<Indicator*>::iterator i = i_w2w1a1.begin(); i != i_w2w1a1.end(); i++)
      output << (*i)->a() << "," << (*i)->b() << "," << (*i)->c() << "," << (*i)->lambda() << endl;
    output.close();

    sst.str("");
    sst << subdir << "/" << "w2w1.csv";
    output.open(sst.str().c_str());
    for(vector<Indicator*>::iterator i = i_w2w1.begin(); i != i_w2w1.end(); i++)
      output << (*i)->a() << "," << (*i)->b() << "," << (*i)->lambda() << endl;
    output.close();

    sst.str("");
    sst << subdir << "/" << "w2a1.csv";
    output.open(sst.str().c_str());
    for(vector<Indicator*>::iterator i = i_w2a1.begin(); i != i_w2a1.end(); i++)
      output << (*i)->a() << "," << (*i)->b() << "," << (*i)->lambda() << endl;
    output.close();

    sst.str("");
    sst << subdir << "/" << "s1w1.csv";
    output.open(sst.str().c_str());
    for(vector<Indicator*>::iterator i = i_s1w1.begin(); i != i_s1w1.end(); i++)
      output << (*i)->a() << "," << (*i)->b() << "," << (*i)->lambda() << endl;
    output.close();

    sst.str("");
    sst << subdir << "/" << "a1s1.csv";
    output.open(sst.str().c_str());
    for(vector<Indicator*>::iterator i = i_a1s1.begin(); i != i_a1s1.end(); i++)
      output << (*i)->a() << "," << (*i)->b() << "," << (*i)->lambda() << endl;
    output.close();
  }

  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int a1 = 0; a1 < bins; a1++)
        m_pw2_c_w1_a1[w2][w1][a1] = pw2_c_w1_a1(w2, w1, a1, i_w2w1a1, i_w2w1, i_w2a1);

  for(int s1 = 0; s1 < bins; s1++)
    for(int a1 = 0; a1 < bins; a1++)
      m_pa1_c_s1[a1][s1] = pa1_c_s1(a1, s1, i_a1s1);

  for(int w1 = 0; w1 < bins; w1++)
    for(int s1 = 0; s1 < bins; s1++)
     m_ps1_c_w1[s1][w1] = ps1_c_w1(s1, w1, i_s1w1);

  for(int w1 = 0; w1 < bins; w1++)
    m_pw1[w1] = pw1(w1);

  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int s1 = 0; s1 < bins; s1++)
        for(int a1 = 0; a1 < bins; a1++)
          pw2w1s1a1[w2][w1][s1][a1] =
            m_pw2_c_w1_a1[w2][w1][a1] *
            m_pa1_c_s1[a1][s1] *
            m_ps1_c_w1[s1][w1] *
            m_pw1[w1];

  if(logLambda)
  {
    stringstream sst;
    sst << "logs/phi_" << phi << "_psi_" << psi << "_chi_" << chi << "_mu_" << mu << "_zeta_" << zeta << "_tau_" << tau;
    string subdir = sst.str();

    sst.str("");
    sst << subdir << "/" << "pw2w1s1a1.csv";
    ofstream output(sst.str().c_str());
    for(int w2 = 0; w2 < bins; w2++)
      for(int w1 = 0; w1 < bins; w1++)
        for(int s1 = 0; s1 < bins; s1++)
          for(int a1 = 0; a1 < bins; a1++)
            output << w2 << "," << w1 << "," << s1 << "," << a1 << "," << pw2w1s1a1[w2][w1][s1][a1] << endl;
    output.close();

    sst.str("");
    sst << subdir << "/" << "pw2_c_w1_a1.csv";
    output.open(sst.str().c_str());
    for(int w2 = 0; w2 < bins; w2++)
      for(int w1 = 0; w1 < bins; w1++)
        for(int a1 = 0; a1 < bins; a1++)
          output << w2 << "," << w1 << "," << a1 << "," << m_pw2_c_w1_a1[w2][w1][a1] << endl;
    output.close();

    sst.str("");
    sst << subdir << "/" << "pa1_c_s1.csv";
    output.open(sst.str().c_str());
    for(int a1 = 0; a1 < bins; a1++)
      for(int s1 = 0; s1 < bins; s1++)
        output << a1 << "," << s1 << "," << m_pa1_c_s1[a1][s1] << endl;
    output.close();

    sst.str("");
    sst << subdir << "/" << "ps1_c_w1.csv";
    output.open(sst.str().c_str());
    for(int s1 = 0; s1 < bins; s1++)
      for(int w1 = 0; w1 < bins; w1++)
        output << s1 << "," << w1 << "," << m_ps1_c_w1[s1][w1] << endl;
    output.close();

    sst.str("");
    sst << subdir << "/" << "pw1.csv";
    output.open(sst.str().c_str());
    for(int w1 = 0; w1 < bins; w1++)
      output << w1 << "," << pw1(w1) << endl;
    output.close();
  }

  i_w2w1a1.clear();
  i_w2w1.clear();
  i_w2a1.clear();
  i_s1w1.clear();
  i_a1s1.clear();
}

double calculate_mc_w(double mu, double phi, double psi, double chi, double zeta, double tau)
{
  double**** m_pw2w1s1a1   = new double***[bins];
  double***  m_pw2_c_w1_a1 = new double**[bins];
  double**   m_pa1_c_s1    = new double*[bins];
  double**   m_ps1_c_w1    = new double*[bins];
  double*    m_pw1         = new double[bins];
  double***  m_pw2w1a1     = new double**[bins];

  create_pw2w1s1a1(m_pw2w1s1a1);
  create_pw2_c_w1a1(m_pw2_c_w1_a1);
  create_pa1_c_s1(m_pa1_c_s1);
  create_ps1_c_w1(m_ps1_c_w1);
  // create_pw1(m_pw1);
  create_pw2w1a1(m_pw2w1a1);

  generate_probability_distrbution(mu, phi, psi, chi, zeta, tau,
                                   m_pw2w1s1a1, m_pw2_c_w1_a1, m_pa1_c_s1, m_ps1_c_w1, m_pw1);


  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int a1 = 0; a1 < bins; a1++)
        m_pw2w1a1[w2][w1][a1] = 0.0;

  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int a1 = 0; a1 < bins; a1++)
        for(int s1 = 0; s1 < bins; s1++)
          m_pw2w1a1[w2][w1][a1] += m_pw2w1s1a1[w2][w1][s1][a1];

  double r = entropy::distribution::MC_W(m_pw2w1a1, bins, bins, bins);

  clear_pw2w1s1a1(m_pw2w1s1a1);
  clear_pw2_c_w1a1(m_pw2_c_w1_a1);
  clear_pa1_c_s1(m_pa1_c_s1);
  clear_ps1_c_w1(m_ps1_c_w1);
  clear_pw1(m_pw1);
  clear_pw2w1a1(m_pw2w1a1);

  return r;
}

double calculate_mc_a(double mu, double phi, double psi, double chi, double zeta, double tau)
{
  double**** m_pw2w1s1a1   = new double***[bins];
  double***  m_pw2_c_w1_a1 = new double**[bins];
  double**   m_pa1_c_s1    = new double*[bins];
  double**   m_ps1_c_w1    = new double*[bins];
  double*    m_pw1         = new double[bins];
  double***  m_pw2a1w1     = new double**[bins];

  create_pw2w1s1a1(m_pw2w1s1a1);
  create_pw2_c_w1a1(m_pw2_c_w1_a1);
  create_pa1_c_s1(m_pa1_c_s1);
  create_ps1_c_w1(m_ps1_c_w1);
  // create_pw1(m_pw1);
  // create_pw2a1w1(m_pw2w1a1);
  create_pw2a1w1(m_pw2a1w1);

  generate_probability_distrbution(mu, phi, psi, chi, zeta, tau,
                                   m_pw2w1s1a1, m_pw2_c_w1_a1, m_pa1_c_s1, m_ps1_c_w1, m_pw1);


  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int a1 = 0; a1 < bins; a1++)
        m_pw2a1w1[w2][a1][w1] = 0.0;

  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int a1 = 0; a1 < bins; a1++)
        for(int s1 = 0; s1 < bins; s1++)
          m_pw2a1w1[w2][a1][w1] += m_pw2w1s1a1[w2][w1][s1][a1];

  double r = entropy::distribution::MC_A(m_pw2a1w1, bins, bins, bins);

  clear_pw2w1s1a1(m_pw2w1s1a1);
  clear_pw2_c_w1a1(m_pw2_c_w1_a1);
  clear_pa1_c_s1(m_pa1_c_s1);
  clear_ps1_c_w1(m_ps1_c_w1);
  clear_pw1(m_pw1);
  clear_pw2a1w1(m_pw2a1w1);

  return r;
}

double calculate_mc_mi(double mu, double phi, double psi, double chi, double zeta, double tau)
{
  double**** m_pw2w1s1a1   = new double***[bins];
  double***  m_pw2_c_w1_a1 = new double**[bins];
  double**   m_pa1_c_s1    = new double*[bins];
  double**   m_ps1_c_w1    = new double*[bins];
  double*    m_pw1         = new double[bins];

  create_pw2w1s1a1(m_pw2w1s1a1);
  create_pw2_c_w1a1(m_pw2_c_w1_a1);
  create_pa1_c_s1(m_pa1_c_s1);
  create_ps1_c_w1(m_ps1_c_w1);

  generate_probability_distrbution(mu, phi, psi, chi, zeta, tau,
                                   m_pw2w1s1a1, m_pw2_c_w1_a1, m_pa1_c_s1, m_ps1_c_w1, m_pw1);

  double r = entropy::distribution::MC_MI(m_pw2w1s1a1, bins, bins, bins, bins);

  clear_pw2w1s1a1(m_pw2w1s1a1);
  clear_pw2_c_w1a1(m_pw2_c_w1_a1);
  clear_pa1_c_s1(m_pa1_c_s1);
  clear_ps1_c_w1(m_ps1_c_w1);
  clear_pw1(m_pw1);

  return r;
}

double calculate_mc_sy_gis(double mu, double phi, double psi, double chi, double zeta, double tau, int iterations)
{
  double**** m_pw2w1s1a1   = new double***[bins];
  double***  m_pw2w1a1     = new double**[bins];
  double***  test_pw2w1a1  = new double**[bins];
  double***  m_pw2_c_w1_a1 = new double**[bins];
  double**   m_pa1_c_s1    = new double*[bins];
  double**   m_ps1_c_w1    = new double*[bins];
  double*    m_pw1         = new double[bins];

  create_pw2w1s1a1(m_pw2w1s1a1);
  create_pw2w1a1(m_pw2w1a1);
  create_pw2w1a1(test_pw2w1a1);
  create_pw2_c_w1a1(m_pw2_c_w1_a1);
  create_pa1_c_s1(m_pa1_c_s1);
  create_ps1_c_w1(m_ps1_c_w1);

  generate_probability_distrbution(mu, phi, psi, chi, zeta, tau,
                                   m_pw2w1s1a1, m_pw2_c_w1_a1, m_pa1_c_s1, m_ps1_c_w1, m_pw1);

  ULContainer* xData = new ULContainer(iterations,2);
  ULContainer* yData = new ULContainer(iterations,1);

  for(int w2 = 0; w2 < bins; w2++)
  {
    for(int w1 = 0; w1 < bins; w1++)
    {
      for(int s1 = 0; s1 < bins; s1++)
      {
        for(int a1 = 0; a1 < bins; a1++)
        {
          m_pw2w1a1[w2][w1][a1] += m_pw2w1s1a1[w2][w1][s1][a1];
        }
      }
    }
  }

  for(int i = 0; i < iterations; i++)
  {
    double p = Random::unit();
    double s = 0.0;
    bool found = false;
    for(int w2 = 0; w2 < bins; w2++)
    {
      for(int w1 = 0; w1 < bins; w1++)
      {
        for(int a1 = 0; a1 < bins; a1++)
        {
          s += m_pw2w1a1[w2][w1][a1];
          if(s >= p)
          {
            (*yData) << w2;
            (*xData) << w1 << a1;
            // cout << "row " << i << ": [" << w1 << "," << a1 << "] -> [" << w2 << "]" << endl;
            test_pw2w1a1[w2][w1][a1] += 1.0;
            found = true;
            break;
          }
          if(found) break;
        }
        if(found) break;
      }
      if(found) break;
    }
  }

  double test = 0.0;
  for(int w2 = 0; w2 < bins; w2++)
  {
    for(int w1 = 0; w1 < bins; w1++)
    {
      for(int a1 = 0; a1 < bins; a1++)
      {
        test_pw2w1a1[w2][w1][a1] /= (float)iterations;
        if(test_pw2w1a1[w2][w1][a1] > 0.0)
        {
          test += m_pw2w1a1[w2][w1][a1] * log2(m_pw2w1a1[w2][w1][a1] / test_pw2w1a1[w2][w1][a1]);
        }
      }
    }
  }

  if(test > 0.1)
  {
    cerr << "Estimation error is: " << test << endl;
  }





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

  double maxError = 10.0;
  while(maxError > 1.0)
  {
    if(independentModel->error() > 1.0) independentModel->iterate();
    if(dependentModel->error()   > 1.0) dependentModel->iterate();
    maxError = MAX(independentModel->error(), dependentModel->error());
  }
  KL* kl = new KL(dependentModel, independentModel);

  double r = kl->divergence2();

  delete kl;
  delete independentModel;
  delete dependentModel;
  delete xData;
  delete yData;

  clear_pw2w1s1a1(m_pw2w1s1a1);
  clear_pw2w1a1(m_pw2w1a1);
  clear_pw2w1a1(test_pw2w1a1);
  clear_pw2_c_w1a1(m_pw2_c_w1_a1);
  clear_pa1_c_s1(m_pa1_c_s1);
  clear_ps1_c_w1(m_ps1_c_w1);
  clear_pw1(m_pw1);

  return r;
}

double calculate_mc_sy_scgis(double mu, double phi, double psi, double chi, double zeta, double tau, int iterations)
{
  double**** m_pw2w1s1a1   = new double***[bins];
  double***  m_pw2w1a1     = new double**[bins];
  double***  test_pw2w1a1  = new double**[bins];
  double***  m_pw2_c_w1_a1 = new double**[bins];
  double**   m_pa1_c_s1    = new double*[bins];
  double**   m_ps1_c_w1    = new double*[bins];
  double*    m_pw1         = new double[bins];

  create_pw2w1s1a1(m_pw2w1s1a1);
  create_pw2w1a1(m_pw2w1a1);
  create_pw2w1a1(test_pw2w1a1);
  create_pw2_c_w1a1(m_pw2_c_w1_a1);
  create_pa1_c_s1(m_pa1_c_s1);
  create_ps1_c_w1(m_ps1_c_w1);

  generate_probability_distrbution(mu, phi, psi, chi, zeta, tau,
                                   m_pw2w1s1a1, m_pw2_c_w1_a1, m_pa1_c_s1, m_ps1_c_w1, m_pw1);

  ULContainer* xData = new ULContainer(iterations,2);
  ULContainer* yData = new ULContainer(iterations,1);

  for(int w2 = 0; w2 < bins; w2++)
  {
    for(int w1 = 0; w1 < bins; w1++)
    {
      for(int s1 = 0; s1 < bins; s1++)
      {
        for(int a1 = 0; a1 < bins; a1++)
        {
          m_pw2w1a1[w2][w1][a1] += m_pw2w1s1a1[w2][w1][s1][a1];
        }
      }
    }
  }

  for(int i = 0; i < iterations; i++)
  {
    double p = Random::unit();
    double s = 0.0;
    bool found = false;
    for(int w2 = 0; w2 < bins; w2++)
    {
      for(int w1 = 0; w1 < bins; w1++)
      {
        for(int a1 = 0; a1 < bins; a1++)
        {
          s += m_pw2w1a1[w2][w1][a1];
          if(s >= p)
          {
            (*yData) << w2;
            (*xData) << w1 << a1;
            test_pw2w1a1[w2][w1][a1] += 1.0;
            found = true;
            break;
          }
          if(found) break;
        }
        if(found) break;
      }
      if(found) break;
    }
  }

  double test = 0.0;
  for(int w2 = 0; w2 < bins; w2++)
  {
    for(int w1 = 0; w1 < bins; w1++)
    {
      for(int a1 = 0; a1 < bins; a1++)
      {
        test_pw2w1a1[w2][w1][a1] /= (float)iterations;
        if(test_pw2w1a1[w2][w1][a1] > 0.0)
        {
          test += m_pw2w1a1[w2][w1][a1] * log2(m_pw2w1a1[w2][w1][a1] / test_pw2w1a1[w2][w1][a1]);
        }
      }
    }
  }

  if(test > 0.1)
  {
    cerr << "Estimation error is: " << test << endl;
  }

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

  double maxError = 10.0;
  while(maxError > 1.0)
  {
    if(independentModel->error() > 1.0) independentModel->iterate();
    if(dependentModel->error()   > 1.0) dependentModel->iterate();
    maxError = MAX(independentModel->error(), dependentModel->error());
  }
  KL* kl = new KL(dependentModel, independentModel);

  double r = kl->divergence2();

  delete kl;
  delete independentModel;
  delete dependentModel;
  delete xData;
  delete yData;

  clear_pa1_c_s1(m_pa1_c_s1);
  clear_ps1_c_w1(m_ps1_c_w1);
  clear_pw1(m_pw1);
  clear_pw2_c_w1a1(m_pw2_c_w1_a1);
  clear_pw2w1a1(m_pw2w1a1);
  clear_pw2w1a1(test_pw2w1a1);
  clear_pw2w1s1a1(m_pw2w1s1a1);

  return r;
}

double calculate_mc_sy_orig_nid(double mu, double phi, double psi, double chi, double zeta, double tau, int iterations)
{
  double**** m_pw2w1s1a1   = new double***[bins];
  double***  m_pw2_c_w1_a1 = new double**[bins];
  double**   m_pa1_c_s1    = new double*[bins];
  double**   m_ps1_c_w1    = new double*[bins];
  double*    m_pw1         = new double[bins];
  double***  m_pw2w1a1     = new double**[bins];

  create_pw2w1s1a1(m_pw2w1s1a1);
  create_pw2_c_w1a1(m_pw2_c_w1_a1);
  create_pa1_c_s1(m_pa1_c_s1);
  create_ps1_c_w1(m_ps1_c_w1);
  create_pw2w1a1(m_pw2w1a1);

  generate_probability_distrbution(mu, phi, psi, chi, zeta, tau,
                                   m_pw2w1s1a1, m_pw2_c_w1_a1, m_pa1_c_s1, m_ps1_c_w1, m_pw1);

  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int a1 = 0; a1 < bins; a1++)
        m_pw2w1a1[w2][w1][a1] = 0.0;

  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int a1 = 0; a1 < bins; a1++)
        for(int s1 = 0; s1 < bins; s1++)
          m_pw2w1a1[w2][w1][a1] += m_pw2w1s1a1[w2][w1][s1][a1];

  vector<vector<int> > features;

  vector<int> a;
  int bits = (int)log2(bins);
  // cout << bits << endl;
  for(int i = 0; i < bits; i++) // W',w
  {
    a.push_back(i);
    a.push_back(i + bits);
  }
  features.push_back(a);

  vector<int> b;
  for(int i = 0; i < bits; i++) // W',A
  {
    b.push_back(i);
    b.push_back(i + 2 * bits);
  }
  features.push_back(b);

  // vector<int> c;
  // for(int i = 0; i < bits; i++) // W,A
  // {
    // c.push_back(i + bits);
    // c.push_back(i + 2 * bits);
  // }
  // features.push_back(c);

  vector<double> p;
  for(int w2 = 0; w2 < bins; w2++)
  {
    for(int w1 = 0; w1 < bins; w1++)
    {
      for(int a1 = 0; a1 < bins; a1++)
      {
        p.push_back(m_pw2w1a1[w2][w1][a1]);
        // cout << w2 << ", " << w1 << ", " << a1 << ": " << m_pw2w1a1[w2][w1][a1] << endl;
      }
    }
  }

  Original* itsc = new Original((int)(log2(p.size())), features, p );
  itsc->iterate(iterations);
  vector<double> pconv = itsc->getp();

  vector<int> x;
  for(int i = 0; i < bits; i++)
  {
    x.push_back(i + bits);
    x.push_back(i + 2 * bits);
  }

  vector<int> y;
  for(int i = 0; i < bits; i++)
  {
    y.push_back(i);
  }

  double r = itsc->calculateConditionalKL(p,pconv,y,x);

  delete itsc;
  // cout << r << endl;
  clear_pw2w1s1a1(m_pw2w1s1a1);
  clear_pw2_c_w1a1(m_pw2_c_w1_a1);
  clear_pa1_c_s1(m_pa1_c_s1);
  clear_ps1_c_w1(m_ps1_c_w1);
  clear_pw1(m_pw1);


  return r;
}

void generate(double mu, double phi, double psi, double chi, double zeta, double tau, int iterations)
{

  double**** m_pw2w1s1a1   = new double***[bins];
  double***  m_pw2w1a1     = new double**[bins];
  double***  test_pw2w1a1  = new double**[bins];
  double***  m_pw2_c_w1_a1 = new double**[bins];
  double**   m_pa1_c_s1    = new double*[bins];
  double**   m_ps1_c_w1    = new double*[bins];
  double*    m_pw1         = new double[bins];

  create_pw2w1s1a1(m_pw2w1s1a1);
  create_pw2w1a1(m_pw2w1a1);
  create_pw2w1a1(test_pw2w1a1);
  create_pw2_c_w1a1(m_pw2_c_w1_a1);
  create_pa1_c_s1(m_pa1_c_s1);
  create_ps1_c_w1(m_ps1_c_w1);

  generate_probability_distrbution(mu, phi, psi, chi, zeta, tau,
                                   m_pw2w1s1a1, m_pw2_c_w1_a1, m_pa1_c_s1, m_ps1_c_w1, m_pw1);

  ULContainer* data = new ULContainer(iterations,3);

  for(int w2 = 0; w2 < bins; w2++)
  {
    for(int w1 = 0; w1 < bins; w1++)
    {
      for(int s1 = 0; s1 < bins; s1++)
      {
        for(int a1 = 0; a1 < bins; a1++)
        {
          m_pw2w1a1[w2][w1][a1] += m_pw2w1s1a1[w2][w1][s1][a1];
        }
      }
    }
  }

  for(int i = 0; i < iterations; i++)
  {
    double p = Random::unit();
    double s = 0.0;
    bool found = false;
    for(int w2 = 0; w2 < bins; w2++)
    {
      for(int w1 = 0; w1 < bins; w1++)
      {
        for(int a1 = 0; a1 < bins; a1++)
        {
          s += m_pw2w1a1[w2][w1][a1];
          if(s >= p)
          {
            (*data) << w2 << w1 << a1;
            test_pw2w1a1[w2][w1][a1] += 1.0;
            found = true;
            break;
          }
          if(found) break;
        }
        if(found) break;
      }
      if(found) break;
    }
  }

  Csv::write("data.csv",data);
  delete data;
}

double calculate_mc_sy_orig(double mu, double phi, double psi, double chi, double zeta, double tau, int iterations)
{
  double**** m_pw2w1s1a1   = new double***[bins];
  double***  m_pw2_c_w1_a1 = new double**[bins];
  double**   m_pa1_c_s1    = new double*[bins];
  double**   m_ps1_c_w1    = new double*[bins];
  double*    m_pw1         = new double[bins];
  double***  m_pw2w1a1     = new double**[bins];

  create_pw2w1s1a1(m_pw2w1s1a1);
  create_pw2_c_w1a1(m_pw2_c_w1_a1);
  create_pa1_c_s1(m_pa1_c_s1);
  create_ps1_c_w1(m_ps1_c_w1);
  create_pw2w1a1(m_pw2w1a1);

  generate_probability_distrbution(mu, phi, psi, chi, zeta, tau,
                                   m_pw2w1s1a1, m_pw2_c_w1_a1, m_pa1_c_s1, m_ps1_c_w1, m_pw1);

  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int a1 = 0; a1 < bins; a1++)
        m_pw2w1a1[w2][w1][a1] = 0.0;

  for(int w2 = 0; w2 < bins; w2++)
    for(int w1 = 0; w1 < bins; w1++)
      for(int a1 = 0; a1 < bins; a1++)
        for(int s1 = 0; s1 < bins; s1++)
          m_pw2w1a1[w2][w1][a1] += m_pw2w1s1a1[w2][w1][s1][a1];

  vector<vector<int> > features;

  vector<int> a;
  int bits = (int)log2(bins);
  // cout << bits << endl;
  for(int i = 0; i < bits; i++) // W',w
  {
    a.push_back(i);
    a.push_back(i + bits);
  }
  features.push_back(a);

  vector<int> b;
  for(int i = 0; i < bits; i++) // W',A
  {
    b.push_back(i);
    b.push_back(i + 2 * bits);
  }
  features.push_back(b);

  vector<int> c;
  for(int i = 0; i < bits; i++) // W,A
  {
    c.push_back(i + bits);
    c.push_back(i + 2 * bits);
  }
  features.push_back(c);

  vector<double> p;
  for(int w2 = 0; w2 < bins; w2++)
  {
    for(int w1 = 0; w1 < bins; w1++)
    {
      for(int a1 = 0; a1 < bins; a1++)
      {
        p.push_back(m_pw2w1a1[w2][w1][a1]);
        // cout << w2 << ", " << w1 << ", " << a1 << ": " << m_pw2w1a1[w2][w1][a1] << endl;
      }
    }
  }

  Original* itsc = new Original((int)(log2(p.size())), features, p );
  itsc->iterate(iterations);
  vector<double> pconv = itsc->getp();

  vector<int> x;
  for(int i = 0; i < bits; i++)
  {
    x.push_back(i + bits);
    x.push_back(i + 2 * bits);
  }

  vector<int> y;
  for(int i = 0; i < bits; i++)
  {
    y.push_back(i);
  }

  double r = itsc->calculateConditionalKL(p,pconv,y,x);

  delete itsc;
  // cout << r << endl;
  clear_pw2w1s1a1(m_pw2w1s1a1);
  clear_pw2_c_w1a1(m_pw2_c_w1_a1);
  clear_pa1_c_s1(m_pa1_c_s1);
  clear_ps1_c_w1(m_ps1_c_w1);
  clear_pw1(m_pw1);


  return r;
}

void deleteAndCreateLogDir()
{
  if(exists(boost::filesystem::path("logs")))
  {
    cout << "removing logs directory." << endl;
    remove_all(boost::filesystem::path("logs"));
  }
  create_directory(boost::filesystem::path("logs"));
}

int main(int argc, char** argv)
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);
  Random::initialise();

  bins        = FLAGS_b;
  logLambda   = FLAGS_L;
  logMatrices = FLAGS_M;

  if(logLambda || logMatrices) deleteAndCreateLogDir();

  if(bins < 2)
  {
    cout << "We need at least two bins per random variable" << endl;
    exit(-1);
  }

  if(FLAGS_o.length() == 0)
  {
    cout << "Please provide an output file." << endl;
    exit(-1);
  }

  cout << "MC:     " << FLAGS_mc   << endl;
  cout << "Output: " << FLAGS_o    << endl;
  cout << "Bins:   " << FLAGS_b    << endl;
  cout << "mu:     " << FLAGS_mu   << endl;
  cout << "phi:    " << FLAGS_phi  << endl;
  cout << "psi:    " << FLAGS_psi  << endl;
  cout << "zeta:   " << FLAGS_zeta << endl;
  cout << "tau:    " << FLAGS_tau  << endl;
  cout << "chi:    " << FLAGS_chi  << endl;

  vector<float> mu   = float_tokenizer(FLAGS_mu);
  vector<float> phi  = float_tokenizer(FLAGS_phi);
  vector<float> psi  = float_tokenizer(FLAGS_psi);
  vector<float> zeta = float_tokenizer(FLAGS_zeta);
  vector<float> tau  = float_tokenizer(FLAGS_tau);
  vector<float> chi  = float_tokenizer(FLAGS_chi);

  // cout << "Phi: " << endl;
  // for(vector<float>::iterator i = phi.begin(); i != phi.end(); i++)
  // {
    // cout << *i << endl;
  // }

  int type = NONE;
  if(FLAGS_mc == "MC_W")        type = __MC_W;
  if(FLAGS_mc == "MC_A")        type = __MC_A;
  if(FLAGS_mc == "MC_MI")       type = __MC_MI;
  if(FLAGS_mc == "MC_SY_GIS")   type = __MC_SY_GIS;
  if(FLAGS_mc == "MC_SY_SCGIS") type = __MC_SY_SCGIS;
  if(FLAGS_mc == "MC_SYO")      type = __MC_SYO;
  if(FLAGS_mc == "MC_SYO_NID")  type = __MC_SYO_NID;
  if(FLAGS_mc == "GENERATE")    type = __GENERATE;
  if(type == NONE)
  {
    cerr << "Error: unknown MC type given \"" << FLAGS_mc << "\"." << endl;
    exit(-1);
  }


  ofstream out;
  if(type != __GENERATE) out.open(FLAGS_o.c_str());

  cout << "Number of iterations: " << mu.size() * phi.size() * psi.size() * zeta.size() * tau.size() * chi.size() << endl;
  boost::progress_display show_progress( phi.size() * psi.size() * chi.size());
  if(type != __GENERATE) out << "# phi, psi, chi, mu, zeta, tau, MC" << endl;

#pragma omp parallel for
  LOOP(i, phi)
  {
    LOOP(j, psi)
    {
      LOOP(n, chi)
      {
        LOOP(k, mu)
        {
          LOOP(l, zeta)
          {
            LOOP(m, tau)
            {
              double r = 0.0;
              switch(type)
              {
                case __MC_W:        r = calculate_mc_w(mu[k],           phi[i], psi[j], chi[n], zeta[l], tau[m]); break;
                case __MC_A:        r = calculate_mc_a(mu[k],           phi[i], psi[j], chi[n], zeta[l], tau[m]); break;
                case __MC_MI:       r = calculate_mc_mi(mu[k],          phi[i], psi[j], chi[n], zeta[l], tau[m]); break;
                case __MC_SY_SCGIS: r = calculate_mc_sy_scgis(mu[k],    phi[i], psi[j], chi[n], zeta[l], tau[m], FLAGS_sysi); break;
                case __MC_SY_GIS:   r = calculate_mc_sy_gis(mu[k],      phi[i], psi[j], chi[n], zeta[l], tau[m], FLAGS_sysi); break;
                case __MC_SYO:      r = calculate_mc_sy_orig(mu[k],     phi[i], psi[j], chi[n], zeta[l], tau[m], FLAGS_syci); break;
                case __MC_SYO_NID:  r = calculate_mc_sy_orig_nid(mu[k], phi[i], psi[j], chi[n], zeta[l], tau[m], FLAGS_syci); break;
                case __GENERATE:    generate(mu[k], phi[i], psi[j], chi[n], zeta[l], tau[m], FLAGS_syci); break;
              }
#pragma omp critical
              {
                if(type != __GENERATE)
                {
                  out
                    << phi[i]  << ","
                    << psi[j]  << ","
                    << chi[n]  << ","
                    << mu[k]   << ","
                    << zeta[l] << ","
                    << tau[m]  << ","
                    << r << endl;
                  out.flush();
                }
              }
            }
          }
        }
#pragma omp critical
              ++show_progress;
      }
    }
  }
  if(type != __GENERATE) out.close();
}
