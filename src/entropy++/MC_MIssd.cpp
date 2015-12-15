#include <entropy++/MC_MIssd.h>

#include <entropy++/MIssd.h>

Container* MC_MI_sparse_matrix_state_dependent(Container* W2, Container* W1, Container* S1, Container* A1, int mode)
{
  Container *mi_w2_w1 = MI_state_dependent_sparse_matrix(W2, W1, mode);
  Container *mi_a1_s1 = MI_state_dependent_sparse_matrix(S1, A1, mode);
  Container *r = new Container(mi_w2_w1->rows(), 1);
  for(int i = 0; i < (int)mi_w2_w1->rows(); i++)
  {
    (*r)(i,0) = mi_w2_w1->get(i,0) - mi_a1_s1->get(i,0);
  }
  return r;
}
