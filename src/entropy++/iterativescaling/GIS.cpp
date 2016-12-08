#include <entropy++/iterativescaling/GIS.h>

using namespace entropy::iterativescaling;

GIS::GIS() : Model()
{ }

void GIS::init()
{
  createUniqueContainer();
  countObservedFeatures();
}

void GIS::iterate()
{
  generateExpected();

  int index = 0;
  cout << "GIS: Nach generateExpected: "<< endl;
  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    cout << "Feature " << index++ << endl;
    int mfindex = 0;
    for(vector<MFeature*>::iterator mf = (*f)->begin(); mf != (*f)->end(); mf++)
    {
      cout << "MF " << mfindex++ << " obs: " << (*mf)->observed() << " exp: " << (*mf)->expected() << " lamda: " << (*mf)->lambda() << endl;
    }
  }

  _error = 0.0;

  double max = 0.0;
  double n   = 0.0;
  double o   = 0.0;


  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    for(vector<MFeature*>::iterator mf = (*f)->begin();
        mf != (*f)->end(); mf++)
    {
      if((*mf)->lambda() > max) max = (*mf)->lambda();
    }
  }

  cout << "Max: " << max << endl;

  cout << "Nach update: "<< endl;
  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    for(vector<MFeature*>::iterator mf = (*f)->begin();
        mf != (*f)->end(); mf++)
    {
      // ueber relations iterieren
      double o = (*mf)->lambda(); // old
      double n = o + 0.1 * (1.0/max)
        * log((*mf)->observed() / (*mf)->expected());
      cout << "o: " << o << " n: " << n << " obs: " << (*mf)->observed() << " exp: " << (*mf)->expected() << endl;
      (*mf)->setLambda(n);
      _error += fabs((*mf)->observed() - (*mf)->expected());
    }
  }
}

double GIS::error()
{
  return _error;
}
