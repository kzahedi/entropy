#include "Random.h"

#include <stdlib.h>

#include <iostream>

#ifdef __linux__
#include <unistd.h>
#include <fcntl.h>
#endif


using namespace std;
using namespace entropy;

void Random::initialise()
{
#ifdef __linux__
  int randomData = open("/dev/random", O_RDONLY);
  int seed = 0;
  read(randomData, &seed, 2);
  cout << "seed: " << seed << endl;
  seed += time(NULL);
  cout << "seed: " << seed << endl;
  srand48(seed);
#else
  sranddev();
#endif

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
