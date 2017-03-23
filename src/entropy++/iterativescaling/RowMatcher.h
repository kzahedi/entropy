#ifndef __ROW_MATCHER_H__
#define __ROW_MATCHER_H__

#include <vector>

using namespace std;


namespace entropy
{
  namespace iterativescaling
  {
    typedef vector<int> Ivector;
    class RowMatcher
    {
      public:
        RowMatcher(int xSize);
        void add_x(int delta_index, int x_row_index);
        void add_y(int delta_index, int y_row_index);

        vector<int>::iterator x_begin(int index);
        vector<int>::iterator x_end(int index);

        vector<int>::iterator y_begin(int index);
        vector<int>::iterator y_end(int index);

      private:
        Ivector* _x_rows;
        Ivector* _y_rows;
        int      _rows;
    };
  }
}

#endif // __ROW_MATCHER_H__
