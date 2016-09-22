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
    Feature(DContainer &aX, DContainer &aY,vector<int>& systX,vector<int>& systY, double valuelambda);
    Feature(DContainer &aX, DContainer &aY,vector<int>& systX,vector<int>& systY, Matrix &lambda);
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
    int    delta(int indexX, int indexY, vector<double> x, vector<double> y);
    double value(vector<double> x,vector<double> y);
    vector<double> index(int index,bool x);
    Feature& operator=(const Feature& c);
    //FeatureMatrix.cpp ITMatrix.cpp IT.cpp FeatureMatrixsp.cpp GIS.cpp GISsp.cpp InstanceMatrix.cpp SCGIS.cpp Test.cpp GISgp.cpp SCGISgp.cpp
  private:

    int 		_sizeDeltaX;
    int			_sizeDeltaY;
    int         _sizeX;
    int         _sizeY;
    vector<int>* _systX;
    vector<int>* _systY;
    DContainer* _X;
    DContainer* _Y;
    Matrix*     _lambda;
};

#endif
