#include <entropy++/iterativescaling/KL.h>

using namespace entropy;
using namespace entropy::iterativescaling;

KL::KL(Model* p, Model *q)
{
  _p = p;
  _q = q;

  cout << "hier 1000" << endl;
  // cout << "calculate probs p" << endl;
  _p->calculateProbabilities();
  cout << "hier 1001" << endl;
  // cout << "calculate probs q" << endl;
  _q->calculateProbabilities();
  cout << "hier 1002" << endl;
  // cout << "done." << endl;

  vector<int> indices_p_x = _p->getAllColumnsForX();
  cout << "hier 1005" << endl;
  vector<int> indices_p_y = _p->getAllColumnsForY();
  cout << "hier 1006" << endl;

  _x_indices = indices_p_x;
  cout << "hier 1007" << endl;
  _y_indices = indices_p_y;
  cout << "hier 1008" << endl;
  //TODO: add _q-indizes
  vector<int> indices_q_x = _q->getAllColumnsForX();
  cout << "hier 1009" << endl;
  vector<int> indices_q_y = _q->getAllColumnsForY();
  cout << "hier 1010" << endl;

  ULContainer* Xalphabet = _p->XAlphabet();
  cout << "hier 1011" << endl;
  ULContainer* Yalphabet = _p->YAlphabet();
  cout << "hier 1012" << endl;
  ULContainer* Xalphabet_both = (Xalphabet->columns(_x_indices))->unique();
  cout << "hier 1013" << endl;
  ULContainer* Yalphabet_both = (Yalphabet->columns(_y_indices))->unique();
  cout << "hier 1014" << endl;
  _xsize = Xalphabet_both->rows();
  cout << "hier 1015" << endl;
  _ysize = Yalphabet_both->rows();
  cout << "hier 1016" << endl;
  ULContainer* x_alphabet_p = (Xalphabet->columns(indices_p_x))->unique();
  cout << "hier 1017" << endl;
  ULContainer* y_alphabet_p = (Yalphabet->columns(indices_p_y))->unique();
  cout << "hier 1018" << endl;
  ULContainer* x_alphabet_q = (Xalphabet->columns(indices_q_x))->unique();
  cout << "hier 1019" << endl;
  ULContainer* y_alphabet_q = (Yalphabet->columns(indices_q_y))->unique();
  cout << "hier 1020" << endl;


  for(int i=0; i< _xsize;i++){
    _x_indices_inAlph_p.push_back(x_alphabet_p->find(Xalphabet_both,i));
    _x_indices_inAlph_q.push_back(x_alphabet_q->find(Xalphabet_both,i));
  }

  for(int i=0; i< _ysize;i++){
    _y_indices_inAlph_p.push_back(y_alphabet_p->find(Yalphabet_both,i));
    _y_indices_inAlph_q.push_back(y_alphabet_q->find(Yalphabet_both,i));
  }

}

double KL::divergence2()
{
  double sum = 0.0;
  double m = 0.0;

  //int nr_of_x   = _x_indices.size(); // Xalphabet
  //int nr_of_y   = _y_indices.size(); // Yalphabet

  for(int x = 0; x < _xsize ; x++)
  {
    for(int y = 0; y < _ysize; y++)
    {
      cout << "hier 000" << endl;
      double p_y_c_x = _p->p_y_c_x(_y_indices_inAlph_p[y],_x_indices_inAlph_p[x]);
      cout << "hier 001" << endl;
      double p_x     = _p->p_x(_x_indices_inAlph_p[x]);
      cout << "hier 002" << endl;
      double q_y_c_x = _q->p_y_c_x(_y_indices_inAlph_q[y],_x_indices_inAlph_q[x]);
      cout << "hier 003" << endl;
      if(p_y_c_x > 0.0 && p_x > 0.0 && q_y_c_x > 0.0)
      {
        m = p_y_c_x * p_x * (log2(p_y_c_x) - log2(q_y_c_x));
        cout << "m: "<< m << endl;
        sum += m;
      }
    }
  }
  return sum;
}

double KL::divergenceN()
{
  double sum = 0.0;

  //int nr_of_x   = _x_indices.size(); // Xalphabet
  //int nr_of_y   = _y_indices.size(); // Yalphabet

  for(int x = 0; x < _xsize ; x++)
  {
    for(int y = 0; y < _ysize; y++)
    {
      double p_y_c_x = _p->p_y_c_x(_y_indices_inAlph_p[y],_x_indices_inAlph_p[x]);
      double p_x     = _p->p_x(_x_indices_inAlph_p[x]);
      double q_y_c_x = _q->p_y_c_x(_y_indices_inAlph_q[y],_x_indices_inAlph_q[x]);
      if(p_y_c_x > 0.0 && p_x > 0.0 && q_y_c_x > 0.0)
      {
        sum += p_y_c_x * p_x * (log(p_y_c_x) - log(q_y_c_x));
      }
    }
  }
  return sum;
}
