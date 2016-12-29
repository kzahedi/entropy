#include "Random.h"

#include <stdlib.h>

#include <iostream>

#ifdef __linux__
#include <linux/random.h>
#endif


using namespace std;
using namespace entropy;

void Random::initialise()
{
#ifdef __linux__
  int randomData = open("/dev/random", O_RDONLY);
  int entropy;
  int result = ioctl(randomData, RNDGETENTCNT, &entropy);
  srand48(result);
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
