#ifndef __FEATURE_H__
#define __FEATURE_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <entropy++/SparseMatrix.h>

#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <math.h>

using namespace std;
using namespace entropy;

namespace entropy
{
  namespace iterativescaling
  {
    class Feature
    {
      public:
        Feature();
        Feature(int sizeX, int sizeY, double valuelambda);
        Feature(int sizeX, int sizeY, SparseMatrix &lambda);
        ~Feature();

        friend std::ostream& operator<<(std::ostream& str,Feature& feature){
          str << "Feature:" << endl;
          for(int i = 0; i < feature._sizeX; i++){
            for(int j = 0; j < feature._sizeY; j++){
              str << feature.getLambda(i,j) << " ";
            }
            str << endl;
          }
          return str;
        };

        int    getLambdaSize();
        double getLambda(int i, int j);
        void   setLambda(int i, int j, double newvalue);
        int    delta(double x,double y,double ax,double ay);
        Feature& operator=(const Feature& c);

      private:
        int           _sizeX;
        int           _sizeY;
        SparseMatrix* _lambda;
    };
  }
}

#endif
