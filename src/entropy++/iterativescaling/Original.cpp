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
  delete _alphabet;
}

void Original::_generateAlphabet(int n){
 // cout << " alphabet " << endl;
  for(int z=0; z< _sizeAlphabet; z++){
    for(int s=n; s> 0; s--){
      (*_alphabet)(z,n-s) = ((int)(z%(int)(pow(2,s))))/(int)(pow(2,s-1));
    }
 //   for(int s=0; s<n;s++){
 //     cout<< (*_alphabet)(z,s);
 //   }
 //   cout << endl;
  }

}
vector<double> Original::getp(){
  return _p1;
}
double Original:: _getprop(vector<double> p, int feat, int ind){
  assert(ind<_sizeAlphabet);
  assert(feat < _features.size());
//  cout << "getprop "<<feat << " " << ind << endl;
  bool found = true;
  double sum = 0.0;
  int featuresize =0;
  for(int i=0; i< _sizeAlphabet; i++){
    featuresize = _features[feat].size();
    for(int j=0; j< featuresize; j++){
      if((*_alphabet)(i,_features[feat][j])!=(*_alphabet)(ind,_features[feat][j])){
        found = false;
      }
    }
    if(found){
      sum+=p[i];
 //     cout<< p[i] <<endl;
    }
    found = true;
  }
//  cout << "sum:  " << sum << endl;
  return sum;
}
double  Original::getMarginalProp(int ind,vector<int> feat, vector<double> p){
  double sum = 0.0;
  bool found = true;
 // cout << " featsize " << feat.size() << endl;
  for(int i=0;i<_sizeAlphabet;i++){
    for(int j=0; j<feat.size();j++){
      if((*_alphabet)(i,feat[j])!=(*_alphabet)(ind,feat[j])){
        found = false;
      }
//      cout << i <<" aufsummieren " << p[i] <<" " << (*_alphabet)(ind,feat[j]) << endl;
    }
    if(found){
      sum+=p[i];
    }
    found = true;
  }
//  cout << "ende marg; " << sum << endl;
  return sum;
}
//p(featMarg|featCond)
double Original::getConditionalProp(vector<int> featMarg, vector<int> featCond, int ind, vector<double> p){
    double sum=0.0;
    featMarg.insert(featMarg.end(), featCond.begin(), featCond.end());
    sum = getMarginalProp(ind, featMarg, p)/getMarginalProp(ind,featCond,p);
   // cout << " Cond: "<< getMarginalProp(ind, featMarg, p) << "/ " << getMarginalProp(ind,featCond,p) << "=" << sum << endl;
    return sum;
}

double Original::calculateKL(int iterations){
  double sum = 0.0;
  if(iterations==0){
    for(int i=0;i< _sizeAlphabet;i++){
      if(_p1[i]>0){
        sum+= _p1[i]*(log2(_p1[i])-log2(_p2[i]));
      }
    }
  }
  else{
    for(int i=0;i< _sizeAlphabet;i++){
      if(_p2[i]>0){
        sum+= _p2[i]*(log2(_p2[i])-log2(_p1[i]));
      }
    }
  }
  return sum;
}
void   Original::iterate(double KL){
  int iterations = 0;
  double kl      = 5;
  while( kl > KL || iterations <=_features.size()){
  //  cout << iterations << " " << kl << endl;
    iterations++;
    int featIndex = (iterations-1)%_features.size();
    if((iterations%2)!=0){
      for(int i=0;i<_sizeAlphabet;i++){
 //        cout <<i << " " << _getprop(_targetp, featIndex,i) << "  " << _p1[i] << " " << _getprop(_p1,featIndex,i) << endl;
        _p2[i]=_getprop(_targetp, featIndex,i)*_p1[i];
        if(_p2[i]!=0){
          _p2[i]=_p2[i]/_getprop(_p1,featIndex,i);
        }
      }
    }
    else{
      for(int i=0;i<_sizeAlphabet;i++){
//        cout <<i << " " << _getprop(_targetp, featIndex,i) << "  " << _p2[i] << " " << _getprop(_p2,featIndex,i) << endl;
        _p1[i]=_getprop(_targetp, featIndex,i)*_p2[i];
        if(_p1[i]!=0){
          _p1[i]=_p1[i]/_getprop(_p2,featIndex,i);
        }
      }
    }
    kl= calculateKL(iterations%2);
//    cout << " iteration: " << it << " p1,p2,p   featindex " << featIndex  << endl;
    for(int i=0; i< _sizeAlphabet;i++){
//      cout << _p1[i] << "  " << _p2[i] << "  " << _targetp[i] << endl;
    }
  }
}
void   Original::iterate(int iterations){
  for(int it=1;it<=iterations;it++){
    int featIndex = (it-1)%_features.size();
    if((it%2)!=0){
      for(int i=0;i<_sizeAlphabet;i++){
 //        cout <<i << " " << _getprop(_targetp, featIndex,i) << "  " << _p1[i] << " " << _getprop(_p1,featIndex,i) << endl;
        _p2[i]=_getprop(_targetp, featIndex,i)*_p1[i];
        if(_p2[i]!=0){
          _p2[i]=_p2[i]/_getprop(_p1,featIndex,i);
        }
      }
    }
    else{
      for(int i=0;i<_sizeAlphabet;i++){
//        cout <<i << " " << _getprop(_targetp, featIndex,i) << "  " << _p2[i] << " " << _getprop(_p2,featIndex,i) << endl;
        _p1[i]=_getprop(_targetp, featIndex,i)*_p2[i];
        if(_p1[i]!=0){
          _p1[i]=_p1[i]/_getprop(_p2,featIndex,i);
        }
      }
    }
//    cout << " iteration: " << it << " p1,p2,p   featindex " << featIndex  << endl;
    for(int i=0; i< _sizeAlphabet;i++){
//      cout << _p1[i] << "  " << _p2[i] << "  " << _targetp[i] << endl;
    }
  }

}

