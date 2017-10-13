#include <gflags/gflags.h>
#include <glog/logging.h>

#include <boost/tokenizer.hpp>

#include <entropy++/distribution/MC_W.h>
#include <entropy++/distribution/MC_A.h>
#include <entropy++/distribution/MC_MI.h>

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

#include <math.h>

#define LOOP(a,b) for(vector<float>::iterator a = b.begin(); a != b.end(); a++)
#define BIN(a) ((a == 0)?-1:1)
#define NONE   -1
#define __MC_W  1001
#define __MC_A  1002
#define __MC_MI 1003

using namespace std;
using namespace boost;
using namespace entropy::distribution;

double**** pw2w1s1a1 = NULL;
double***  pw2w1a1   = NULL;
double***  pw2a1w1   = NULL;

vector<float> float_tokenizer(string s)
{
  boost::char_separator<char> sep(":");
  tokenizer<boost::char_separator<char> > tk(s,sep);
  vector<float> vec;
  for(tokenizer<boost::char_separator<char> >::iterator i(tk.begin()); i!=tk.end();++i)
  {
    vec.push_back(atof((*i).c_str()));
  }

  double start = vec[0];
  double delta = vec[1];
  double end   = vec[2];

  double value = start;

  vector<float> values;
  while(value <= end)
  {
    values.push_back(value);
    value += delta;
  }

  return values;
}

void clear_pw2w1s1a1()
{
  if(pw2w1s1a1 != NULL)
  {
    for(int i = 0; i < 2; i++)
    {
      for(int j = 0; j < 2; j++)
      {
        for(int k = 0; k < 2; k++)
        {
          delete pw2w1s1a1[i][j][k];
        }
        delete pw2w1s1a1[i][j];
      }
      delete pw2w1s1a1[i];
    }
    delete pw2w1s1a1;
  }
  pw2w1s1a1 = NULL;
}

void create_pw2w1s1a1()
{
  pw2w1s1a1 = new double***[2];
  for(int i = 0; i < 2; i++)
  {
    pw2w1s1a1[i] = new double**[2];
    for(int j = 0; j < 2; j++)
    {
      pw2w1s1a1[i][j] = new double*[2];
      for(int k = 0; k < 2; k++)
      {
        pw2w1s1a1[i][j][k] = new double[2];
      }
    }
  }
}

void create_pw2w1a1()
{
  pw2w1a1 = new double**[2];
  for(int i = 0; i < 2; i++)
  {
    pw2w1a1[i] = new double*[2];
    for(int j = 0; j < 2; j++)
    {
      pw2w1a1[i][j] = new double[2];
    }
  }
}

void clear_pw2w1a1()
{
  if(pw2w1a1 != NULL)
  {
    for(int i = 0; i < 2; i++)
    {
      for(int j = 0; j < 2; j++)
      {
        delete pw2w1a1[i][j];
      }
      delete pw2w1a1[i];
    }
    delete pw2w1a1;
  }
  pw2w1a1 = NULL;
}

void create_pw2a1w1()
{
  pw2a1w1 = new double**[2];
  for(int i = 0; i < 2; i++)
  {
    pw2a1w1[i] = new double*[2];
    for(int j = 0; j < 2; j++)
    {
      pw2a1w1[i][j] = new double[2];
    }
  }
}

void clear_pw2a1w1()
{
  if(pw2a1w1 != NULL)
  {
    for(int i = 0; i < 2; i++)
    {
      for(int j = 0; j < 2; j++)
      {
        delete pw2a1w1[i][j];
      }
      delete pw2a1w1[i];
    }
    delete pw2a1w1;
  }
  pw2a1w1 = NULL;
}


double pw2_c_w1_a1(int w2, int w1, int a1, double phi, double psi)
{
  double z = exp(phi * BIN(w2) * BIN(w1) + psi * BIN(w2) * BIN(a1));
  double n = exp(phi * -1 * BIN(w1) + psi * -1 * BIN(a1)) +
             exp(phi *  1 * BIN(w1) + psi *  1 * BIN(a1));
  return z/n;
}

double pa1_c_s1(int a1, int s1, double mu)
{
  double z = exp(mu * BIN(a1) * BIN(s1));
  double n = exp(mu * BIN(s1)) + exp(-mu * BIN(s1));
  return z/n;
}

double ps1_c_w1(int s1, int w1, double zeta)
{
  double z = exp(zeta * BIN(w1) * BIN(s1));
  double n = exp(zeta * BIN(w1)) + exp(-zeta * BIN(w1));
  return z/n;
}

double pw1(int w1, double tau)
{
  double z = exp(tau * BIN(w1));
  double n = exp(tau) + exp(-tau);
  return z/n;
}

string toString()
{
  stringstream sst;
  sst << "# w2, w1, s1, a1, value" << endl;
  for(int w2 = 0; w2 < 2; w2++)
    for(int w1 = 0; w1 < 2; w1++)
      for(int s1 = 0; s1 < 2; s1++)
        for(int a1 = 0; a1 < 2; a1++)
          sst << BIN(w2) << "," << BIN(w1) << "," << BIN(s1) << "," << BIN(a1) << ","
            << pw2w1s1a1[w2][w1][s1][a1] << endl;
  return sst.str();
};

void generate_probability_distrbution_pw2w1a1(double mu, double phi, double psi, double zeta, double tau)
{
  // clear_pw2w1s1a1();

  for(int w2 = 0; w2 < 2; w2++)
    for(int w1 = 0; w1 < 2; w1++)
      for(int s1 = 0; s1 < 2; s1++)
        for(int a1 = 0; a1 < 2; a1++)
          pw2w1s1a1[w2][w1][s1][a1] =
            pw2_c_w1_a1(w2, w1, a1, phi, psi) *
            pa1_c_s1(a1, s1, mu) *
            ps1_c_w1(s1, w1, zeta) *
            pw1(w1, tau);
}

DEFINE_string(mu,   "0", "mu values");
DEFINE_string(phi,  "0.0:0.1:5.0", "phi values");
DEFINE_string(psi,  "0.0:0.1:5.0", "psi values");
DEFINE_string(zeta, "0.0", "zeta values");
DEFINE_string(tau,  "0.0", "tau values");
DEFINE_string(mc,   "MC_W",    "quantification");
DEFINE_string(o,    "out.csv",        "output file");

double calculate_mc_w()
{
  for(int w2 = 0; w2 < 2; w2++)
    for(int w1 = 0; w1 < 2; w1++)
      for(int a1 = 0; a1 < 2; a1++)
        pw2w1a1[w2][w1][a1] = 0.0;

  for(int w2 = 0; w2 < 2; w2++)
    for(int w1 = 0; w1 < 2; w1++)
      for(int a1 = 0; a1 < 2; a1++)
        for(int s1 = 0; s1 < 2; s1++)
          pw2w1a1[w2][w1][a1] += pw2w1s1a1[w2][w1][s1][a1];

  return entropy::distribution::MC_W(pw2w1a1, 2, 2, 2);
}

double calculate_mc_a()
{
  for(int w2 = 0; w2 < 2; w2++)
    for(int w1 = 0; w1 < 2; w1++)
      for(int a1 = 0; a1 < 2; a1++)
        pw2a1w1[w2][a1][w1] = 0.0;

  for(int w2 = 0; w2 < 2; w2++)
    for(int w1 = 0; w1 < 2; w1++)
      for(int a1 = 0; a1 < 2; a1++)
        for(int s1 = 0; s1 < 2; s1++)
          pw2a1w1[w2][a1][w1] += pw2w1s1a1[w2][w1][s1][a1];

  return entropy::distribution::MC_A(pw2a1w1, 2, 2, 2);
}

double calculate_mc_mi()
{
  return entropy::distribution::MC_MI(pw2w1s1a1, 2, 2, 2, 2);
}

int main(int argc, char** argv)
{
  google::ParseCommandLineFlags(&argc, &argv, true);
  google::InitGoogleLogging(argv[0]);

  if(FLAGS_o.length() == 0)
  {
    cout << "Please provide an output file." << endl;
    exit(-1);
  }

  cout << "mu:   " << FLAGS_mu   << endl;
  cout << "phi:  " << FLAGS_phi  << endl;
  cout << "psi:  " << FLAGS_psi  << endl;
  cout << "zeta: " << FLAGS_zeta << endl;
  cout << "tau:  " << FLAGS_tau  << endl;

  vector<float> mu   = float_tokenizer(FLAGS_mu);
  vector<float> phi  = float_tokenizer(FLAGS_phi);
  vector<float> psi  = float_tokenizer(FLAGS_psi);
  vector<float> zeta = float_tokenizer(FLAGS_zeta);
  vector<float> tau  = float_tokenizer(FLAGS_tau);

  int type = NONE;
  if(FLAGS_mc == "MC_W") type = __MC_W;
  if(FLAGS_mc == "MC_A") type = __MC_A;
  if(FLAGS_mc == "MC_MI") type = __MC_MI;

  create_pw2w1s1a1();
  switch(type)
  {
    case __MC_W: create_pw2w1a1(); break;
    case __MC_A: create_pw2a1w1(); break;
  }

  ofstream out(FLAGS_o.c_str());

  LOOP(i_mu, mu)
  {
    LOOP(i_phi, phi)
    {
      LOOP(i_psi, psi)
      {
        LOOP(i_zeta, zeta)
        {
          LOOP(i_tau, tau)
          {
            generate_probability_distrbution_pw2w1a1(*i_mu, *i_phi, *i_psi, *i_zeta, *i_tau);
            double r = 0.0;
            switch(type)
            {
              case __MC_W:  r = calculate_mc_w();  break;
              case __MC_A:  r = calculate_mc_a();  break;
              case __MC_MI: r = calculate_mc_mi(); break;
            }
            out << *i_mu << "," << *i_phi << "," << *i_psi << "," << *i_zeta << "," << *i_tau << "," << r << endl;
            out.flush();
          }
        }
      }
    }
  }
  out.close();
  clear_pw2w1s1a1();
  switch(type)
  {
    case __MC_W: clear_pw2w1a1(); break;
    case __MC_A: clear_pw2a1w1(); break;
  }
}
