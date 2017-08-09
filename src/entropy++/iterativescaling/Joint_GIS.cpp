#include <entropy++/defs.h>
#include <entropy++/iterativescaling/Joint_GIS.h>

//#include <glog/logging.h>

#ifdef USE_OPENMP
#  include <omp.h>
#endif

// #include <glog/logging.h>

#define MIN_S 0.0000001

using namespace entropy::iterativescaling;

Joint_GIS::Joint_GIS() : IS()
{}

Joint_GIS::~Joint_GIS()
{
  delete _deltaMatcher;
}
void Joint_GIS::setFeatures(vector<vector<int> > Indices,  vector<Feature*> f)
{
  setXIndices(Indices);
  features = f;
  _Indices=Indices;
}

void Joint_GIS::countObservedFeatures()
{
  for (int i = 0; i < Xdata->rows(); i++)
  {
    vector<unsigned long> xrow = Xdata->row(i);
    for (vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
    {
      bool found = false;
      for (vector<Delta*>::iterator d = (*f)->begin(); d != (*f)->end(); d++)
      {
        if ((*d)->matchX(xrow))
        {
          (*d)->incObserved();
          found = true;
        }
      }
      if (found == false)
      {
        int xI = (*f)->xListIndex();
        vector<int> xIndices;
        vector<int> yIndices;
        if(xI >= 0) xIndices = _Indices[(*f)->xListIndex()];
        vector<unsigned long> yrow;
        Delta *d = new Delta(xrow,xIndices, yrow, yIndices);
        d->setInputOnly();
        deltas.push_back(d);
        (*f)->push_back(d);
      }
    }
  }
  // normalizing observed:
  for (vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    for (vector<Delta*>::iterator d = (*f)->begin(); d != (*f)->end(); d++)
    {
      double c =((*d)->observed())/Xdata->rows();
      (*d)->setObserved(c);
    }
  }
}
void Joint_GIS::setAlphabet(ULContainer* X){
  Xalphabet = X;
}
void Joint_GIS::init()
{
  countObservedFeatures();

  // cout << "Initialising DeltaMatcher" << endl;
  _deltaMatcher = new DeltaMatcher(Xalphabet->rows());
  // boost::progress_display show_progress( Xdata->rows() );
  for(int x = 0; x < Xalphabet->rows(); x++)
  {
    vector<unsigned long> x_row = Xalphabet->row(x);
    for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchX(x_row))
        {
          bool found = false;
          for(vector<Delta*>::iterator dm = _deltaMatcher->begin(x); dm != _deltaMatcher->end(x); dm++)
          {
            if(*dm == *d)
            {
              found = true;
              break;
            }
          }
          if(found == false) _deltaMatcher->add(x, *d);
        }
      }
    }
    // ++show_progress;
}

void Joint_GIS::iterate()
{
  for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  {
    (*d)->setExpected(0.0);
  }
  double z = 0.0;
  vector<double> s;
  for(int j = 0; j < Xalphabet->rows(); j++)
  {
    vector<unsigned long> x_row = Xalphabet->row(j);
    s.resize(Xalphabet->rows());
    for(vector<Delta*>::iterator d = _deltaMatcher->begin(j); d != _deltaMatcher->end(j); d++)
    {
      if((*d)->matchX(x_row))
      {
        s[j] += (*d)->lambda();
      }
    }
    z+= exp(s[j]);
  }
  for(int j = 0; j < Xalphabet->rows(); j++)
  {
    vector<unsigned long> x_row = Xalphabet->row(j);
    for(vector<Delta*>::iterator d = _deltaMatcher->begin(j); d != _deltaMatcher->end(j); d++)
    {
      if((*d)->matchX(x_row))
      {
        (*d)->setExpected((*d)->expected() + exp(s[j]) / z);
      }
    }
  }
  _error = 0.0;

  double max = 1.0;
  double e = 0.0;

  // get f^#
  for(int j = 0; j < Xalphabet->rows(); j++)
  {
    vector<unsigned long> x_row = Xalphabet->row(j);
    {
      double f = 0.0;
      for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchX(x_row))
        {
// #pragma omp atomic update
          f += 1.0;
        }
      }
      if(f > max) max = f;
    }
  }

  for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  {
    double n = (*d)->lambda() + fabs(1.0/max) * log((*d)->observed() / (*d)->expected());
    (*d)->setLambda(n);
    e = (*d)->observed() - (*d)->expected();
    _error += e * e;
  }

  _error = sqrt(_error);

}
void Joint_GIS::calculateProbabilities()
{
  vector<double> M(Xalphabet->rows());

  for(int x = 0; x < Xalphabet->rows(); x++)
  {
    vector<unsigned long> x_row = Xalphabet->row(x);
    for (vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
    {
       if((*d)->matchX(x_row))
        {
          M[x] += (*d)->lambda();
        }
     }
   }

  for (int x = 0; x < Xalphabet->rows(); x++)
  {
    M[x]= exp(M[x]);
  }
  double z = 0;
  for (int x = 0; x < Xalphabet->rows(); x++)
  {
    z+= M[x];
  }
  _marginals.resize(Xalphabet->rows());
  for (int x = 0; x < Xalphabet->rows(); x++)
  {
    _marginals[x]= M[x]/z;
  }
  for(int i=0;i< Xalphabet->rows();i++){
    cout<< _marginals[i] << " " << endl;
  }
}
double Joint_GIS::error()
{
  return _error;
}
