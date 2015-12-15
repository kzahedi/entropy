#include <entropy++/MC_Ws.h>

#include <entropy++/sparse/CMI.h>

using namespace entropy::sparse;

double MC_W_sparse_matrix(Container* W2, Container* W1, Container* A1, int mode)
{
  return CMI(W2, W1, A1, mode);
}
