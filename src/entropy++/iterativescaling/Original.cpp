#include "Original.h"
#include <entropy++/powi.h>

#include <iostream>

using namespace entropy::iterativescaling;
using namespace std;

Original::Original()
{
  _sizeAlphabet = 0;
  _alphabet     = NULL;
}

Original::Original(int n, vector<vector<int> > features, vector<double> p)
{
  cout << "Konstruktor" << endl;
  _sizeAlphabet = pow(2,n);
  _alphabet = new Matrix(_sizeAlphabet,n);
  _generateAlphabet(n);

  assert(p.size()==_sizeAlphabet);
  _targetp  = p;

  _features = features;

  _p1 = vector<double>(_sizeAlphabet);
  _p2 = vector<double>(_sizeAlphabet);

  double val = 1.0/_sizeAlphabet;
  for(int i=0; i<_sizeAlphabet;i++){
    _p1[i] = val;
  }
}

Original::~Original()
{
}

void Original::_generateAlphabet(int n){
  cout << " alphabet " << endl;
  for(int z=0; z< _sizeAlphabet; z++){
    for(int s=n; s> 0; s--){
      (*_alphabet)(z,s-1) = ((int)(z%(int)(pow(2,s))))/(int)(pow(2,s-1));
       cout << (*_alphabet)(z,s-1);
    }
    cout << endl;
  }
}
vector<double> Original::getp(){
  return _p1;
}
double Original:: _getprop(vector<double> p, int feat, int ind){
  assert(ind<_sizeAlphabet);
  assert(feat < _features.size());
  double sum = 0.0;
  int featuresize =0;
  for(int i=0; i< _sizeAlphabet; i++){
    featuresize = _features[feat].size();
    for(int j=0; j< featuresize; j++){
      if((*_alphabet)(i,_features[feat][j])==(*_alphabet)(ind,_features[feat][j])){
        sum+=p[i];
      }
    }
  }
  return sum;
}
double  Original::getMarginalProp(int ind,int feat, vector<double> p){
  assert(feat<log2(_sizeAlphabet));
  double sum = 0.0;
  for(int i=0;i<_sizeAlphabet;i++){
    if((*_alphabet)(i,feat)==(*_alphabet)(ind,feat)){
//      cout << i <<" aufsummieren " << p[i] <<" " << (*_alphabet)(ind,feat) << endl;
      sum+=p[i];
    }
  }
//  cout << "ende marg" << endl;
  return sum;
}
void   Original::iterate(int iterations){
  for(int it=1;it<=iterations;it++){
    int featIndex = (it-1)%_features.size();
    if((it%2)!=0){
      for(int i=0;i<_sizeAlphabet;i++){
    //    cout <<i << " " << _getprop(_targetp, featIndex,i) << "  " << _p1[i] << " " << _getprop(_p1,featIndex,i) << endl;
        _p2[i]=_getprop(_targetp, featIndex,i)*_p1[i];
        if(_p2[i]!=0){
          _p2[i]=_p2[i]/_getprop(_p1,featIndex,i);
        }
      }
    }
    else{
      for(int i=0;i<_sizeAlphabet;i++){
  //      cout <<i << " " << _getprop(_targetp, featIndex,i) << "  " << _p2[i] << " " << _getprop(_p2,featIndex,i) << endl;
        _p1[i]=_getprop(_targetp, featIndex,i)*_p2[i];
        if(_p1[i]!=0){
          _p1[i]=_p1[i]/_getprop(_p2,featIndex,i);
        }
      }
    }
 //   cout << " iteration: " << it << " p1,p2,p   featindex " << featIndex  << endl;
    for(int i=0; i< _sizeAlphabet;i++){
 //     cout << _p1[i] << "  " << _p2[i] << "  " << _targetp[i] << endl;
    }
  }

}

