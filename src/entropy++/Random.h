#ifndef __RANDOM_H__
#define __RANDOM_H__

namespace entropy
{
  class Random
  {
    public:
      static void   initialise();
      static void   initialise(int seed);
      static double unit();
      static int    rand(int min, int max);
  };
}

#endif // __RANDOM_H__
