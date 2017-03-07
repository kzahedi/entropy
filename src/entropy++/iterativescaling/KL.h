#ifndef __KL_H__
#define __KL_H__

#include <entropy++/iterativescaling/Model.h>


namespace entropy
{
  namespace iterativescaling
  {
    class KL
    {
      public:
        KL(Model* p, Model *q);

        double divergence2();
        double divergenceN();

      private:
        Model* _p;
        Model* _q;
        vector<int> _x_indices;
        vector<int> _y_indices;
        vector<int> _x_indices_inAlph_p;
        vector<int> _y_indices_inAlph_p;
        vector<int> _x_indices_inAlph_q;
        vector<int> _y_indices_inAlph_q;
        int         _xsize;
        int         _ysize;
    };
  }
}

#endif // __KL_H__
