#ifndef __MIssd_H__
#define __MIssd_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

Container* MI_state_dependent_sparse_matrix(Container* X, Container* Y, int mode = EMPERICAL);

# define MIssd(x,y)   MI_state_dependent_sparse_matrix(x,y)

#endif // __MIssd_H__
