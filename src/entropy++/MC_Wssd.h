#ifndef __MC_Wssd_H__
#define __MC_Wssd_H__

#include <entropy++/Entropy.h>

Container* MC_W_sparse_matrix_state_dependent(Container* W2,
                                              Container* W1,
                                              Container* A1,
                                              int mode = EMPERICAL);

#define MC_Wssd(x, y, z) MC_W_sparse_matrix_state_dependent(x, y, z)

#endif // __MC_W_H__
