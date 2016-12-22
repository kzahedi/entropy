#include "Random.h"

#include <stdlib.h>

#include <iostream>

using namespace std;
using namespace entropy;

void Random::initialise()
{
  time_t t;
  time(&t);
  srand48(t);

  cout << "random initialised:";
  for(int i = 0; i < 10; i++)
  {
    cout << " " << rand(0, 100);
  }
  cout << endl;
  cout << "random initialised:";
  for(int i = 0; i < 10; i++)
  {
    cout << " " << unit();
  }
  cout << endl;
}

double Random::unit()
{
  return drand48();
}

void Random::initialise(int seed)
{
  cout << "using Random::initialise(" << seed << ")" << endl;
  srand(seed);
}

int Random::rand(int min, int max)
{
  return min + int(drand48() * (double)max + 0.5);
}
