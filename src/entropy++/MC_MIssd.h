#ifndef __MC_MI_STATE_DEPENDENT_SPARSE_MATRIX_H__
#define __MC_MI_STATE_DEPENDENT_SPARSE_MATRIX_H__

#include <entropy++/Container.h>

Container* MC_MI_sparse_matrix_state_dependent(Container* W2, Container* W1, Container* S1, Container* A1, int mode = EMPERICAL);

#define MC_MIssd(a,b,c,d) MC_MI_sparse_matrix_state_dependent(a,b,c,d)

#endif
