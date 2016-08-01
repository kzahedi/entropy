#ifndef __DEFS_H__
#define __DEFS_H__

# define EMPERICAL 2001
# define UNIFORM   1001

#define ASSERT(condition, message) \
    if (! (condition)) { \
      std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
      << " line " << __LINE__ << ": " << message << std::endl; \
      std::exit(EXIT_FAILURE); \
    }

#endif // __DEFS_H__
