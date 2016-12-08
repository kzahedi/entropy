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

  _error = 0.0;

  double max = 0.0;
  double n   = 0.0;
  double o   = 0.0;

  for(vector<MFeature*>::iterator f = features.begin();
      f != features.end(); f++)
  {
    if((*f)->lambda() > max) max = (*f)->lambda();
  }

  // ueber relations iterieren
  for(vector<MFeature*>::iterator f = features.begin();
      f != features.end(); f++)
  {
    double o = (*f)->lambda(); // old
    double n = o + 0.1 * (1.0/max)
      * log((*f)->observed() / (*f)->expected());
    // cout << o << " -> " << n << endl;
    (*f)->setLambda(n);
    _error += fabs((*f)->observed() - (*f)->expected());
  }

  // cout << "1: " << features[0]->observed() << " - " << features[0]->expected() << " -> " << features[0]->observed() - features[0]->expected() << endl;


  // cout << "2: " << features[0]->observed() << " - " << features[0]->expected() << " -> " << features[0]->observed() - features[0]->expected() << endl;

}

double GIS::error()
{
  return _error;
}
