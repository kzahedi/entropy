#ifndef __FEATURE_H__
#define __FEATURE_H__

#include <entropy++/Container.h>
#include <entropy++/defs.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <math.h>
#include <entropy++/Matrix.h>

using namespace std;

class Feature
{
  public:
    Feature();
    Feature(DContainer &aX, DContainer &aY,int colValX, int colValY, int systXsize,int systYsize , double valuelambda);
    Feature(DContainer &aX, DContainer &aY,int colValX, int colValY, int sizeSystX,int sizeSystY, Matrix &lambda);
    ~Feature();

    friend std::ostream& operator<<(std::ostream& str,Feature& feature){
      str<< "Feature:" <<endl;
      for(int i=0; i<feature._sizeDeltaX; i++){
        for(int j=0; j<feature._sizeDeltaY; j++){
          str<< feature.getlambda(i,j) << " ";
        }
        str<< endl;
      }
      return str;
    };

    double getlambda(int i, int j);
    void   setlambda(int i, int j, double newvalue);
    int    delta(double x,double y,double ax,double ay);
    Feature& operator=(const Feature& c);
    //evtl. return index

  private:
    int 		_sizeDeltaX;
    int			_sizeDeltaY;
    int         _sizeX;
    int         _sizeY;
    DContainer* _X;
    DContainer* _Y;
    Matrix*     _lambda;
};

#endif
