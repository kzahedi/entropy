#include <entropy++/MC_Wssd.h>

#include <entropy++/CMIssd.h>

Container* MC_W_sparse_matrix_state_dependent(Container* W2, Container* W1, Container* A1, int mode)
{
  return CMI_sparse_matrix_state_dependent(W2, W1, A1, mode);
}
