#include <entropy++/iterativescaling/KL.h>

using namespace entropy;
using namespace entropy::iterativescaling;

KL::KL(Model* p, Model *q)
{
  _p = p;
  _q = q;

  _p->calculateProbabilities();
  _q->calculateProbabilities();
}

double KL::divergence2()
{
  double sum = 0.0;

  int nr_of_x = _p->getNrOfUniqueX(); // Xalphabet
  int nr_of_y = _p->getNrOfUniqueY(); // Yalphabet

#pragma omp parallel for reduction(+:sum)
  for(int x = 0; x < nr_of_x; x++)
  {
    for(int y = 0; y < nr_of_y; y++)
    {
      double p_y_c_x = _p->p_y_c_x_d(y,x);
      double p_x     = _p->p_x_d(x);
      double q_y_c_x = _q->p_y_c_x_d(y,x);
      if(p_y_c_x > 0.0 && p_x > 0.0 && q_y_c_x > 0.0)
      {
        sum += p_y_c_x * p_x * (log2(p_y_c_x) - log2(q_y_c_x));
      }
    }
  }
  return sum;
}

double KL::divergenceN()
{
  double sum = 0.0;

  int nr_of_x = _p->getNrOfUniqueX();
  int nr_of_y = _p->getNrOfUniqueY();

#pragma omp parallel for reduction(+:sum)
  for(int x = 0; x < nr_of_x; x++)
  {
    for(int y = 0; y < nr_of_y; y++)
    {
      double p_y_c_x = _p->p_y_c_x_d(y,x);
      double p_x     = _p->p_x_d(x);
      double q_y_c_x = _q->p_y_c_x_d(y,x);
      if(p_y_c_x > 0.0 && p_x > 0.0 && q_y_c_x > 0.0)
      {
        sum += p_y_c_x * p_x * (log(p_y_c_x) - log(q_y_c_x));
      }
    }
  }
  return sum;
}
