#ifndef __DEFS_H__
#define __DEFS_H__

#include <vector>
#include <string>
#include <iostream>
using namespace std;

#define MAX(a,b)    ((a>b)?a:b)
#define MAX3(a,b,c) MAX(a, MAX(b,c))

#define MIN(a,b)    ((a<b)?a:b)
#define MIN3(a,b,c) MIN(a,MIN(b,c))


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

#ifndef ROW_OUTPUT
#define ROW_OUTPUT(label, vec) \
  cout << label << ":"; \
  for(vector<unsigned long>::iterator i = vec.begin(); i != vec.end(); i++)\
  {\
    cout << " " << *i;\
  }\
  cout << endl;
#endif

#endif // __DEFS_H__
