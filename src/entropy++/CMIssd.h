#ifndef __CMIssd_H__
#define __CMIssd_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>

Container* CMI_sparse_matrix_state_dependent(Container* X, Container* Y, Container* Z, int mode = EMPERICAL);

#define CMIssd(x,y,z) CMI_sparse_matrix_state_dependent(x,y,z)

#endif // __CMIssd_H__
