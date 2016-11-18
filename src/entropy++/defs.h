#ifndef __DEFS_H__
#define __DEFS_H__

#include <vector>
using namespace std;

# define EMPERICAL 2001
# define UNIFORM   1001

#define MC_CW_MODE_ENTROPY 1001
#define MC_CW_MODE_MI      1002


#define ASSERT(condition, message) \
    if (! (condition)) { \
      std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
      << " line " << __LINE__ << ": " << message << std::endl; \
      std::exit(EXIT_FAILURE); \
    }

typedef vector<int>          ivector;
typedef vector<vector<int> > ivvector;

typedef vector<double>          dvector;
typedef vector<vector<double> > dvvector;



#endif // __DEFS_H__

