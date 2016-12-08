#include <entropy++/iterativescaling/Model.h>

using namespace entropy;
using namespace entropy::iterativescaling;

Model::Model()
{
}

void Model::setData(ULContainer *X, ULContainer *Y)
{
  _nrX = X->columns();
  _nrY = Y->columns();

  _X = X;
  _Y = Y;

}

Model::~Model()
{
  delete _X;
  delete _Y;

  for(int i = 0; i < _Xindices.size(); i++)
  {
    delete _uniqueX[i];
  }

  for(int i = 0; i < _Yindices.size(); i++)
  {
    delete _uniqueY[i];
  }

  delete _uniqueX;
  delete _uniqueY;
}



void Model::setFeatures(vector<vector<int> > Xindices,
                        vector<vector<int> > Yindices,
                        vector<Feature*>     f)
{
  _Xindices = Xindices;
  _Yindices = Yindices;
  features  = f;
}

void Model::createUniqueContainer()
{
  _uniqueX = new ULContainer*[_Xindices.size()];
  _uniqueY = new ULContainer*[_Yindices.size()];

  for(int i = 0; i < _Xindices.size(); i++)
  {
    _uniqueX[i] = _X->columns(_Xindices[i]);
    _uniqueX[i] = _uniqueX[i]->unique();
  }
  for(int i = 0; i < _Yindices.size(); i++)
  {
    _uniqueY[i] = _Y->columns(_Yindices[i]);
    _uniqueY[i] = _uniqueY[i]->unique();
  }
}

void Model::countObservedFeatures()
{
  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    int i = (*f)->xListIndex();
    int k = (*f)->yListIndex();
    ULContainer* xMatchColumns = _X->columns(_Xindices[i]);
    ULContainer* yMatchColumns = _Y->columns(_Yindices[k]);
    for(int d = 0; d < _X->rows(); d++)
    {
      int j = _uniqueX[i]->find(xMatchColumns,d);
      int l = _uniqueY[k]->find(yMatchColumns,d);

      bool found = false;
      for(vector<Delta*>::iterator mf = (*f)->begin(); mf != (*f)->end(); mf++)
      {
        if((*mf)->match(j,l))
        {
          (*mf)->incObserved();
          found = true;
          break;
        }
      }
      if(found == false)
      {
        Delta *mf = new Delta(j,l);
        mf->incObserved();
        (*f)->push_back(mf);
      }
    }
    delete xMatchColumns;
    delete yMatchColumns;
  }

  // obs(x)
  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    int i = (*f)->xListIndex(); // list of columns that define X
    ULContainer* xUniqueForThisFeature = _uniqueX[i];
    ULContainer* xMatchColumns         = _X->columns(_Xindices[i]);
    for(int x = 0; x < xUniqueForThisFeature->rows(); x++)
    {
      int count = xMatchColumns->findlist(xUniqueForThisFeature,x).size();
      (*f)->setUniqueXCount(x, count);
    }
  }

  // int index = 0;
  // cout << "In countObservedFeatures: "<< endl;
  // for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  // {
    // cout << "Feature " << index++ << endl;
    // int i = (*f)->xListIndex(); // list of columns that define X
    // for(int x = 0; x < _uniqueX[i]->rows(); x++)
    // {
      // cout << "  " << x << " x count: " << (*f)->getUniqueXCount(x) << endl;
    // }
    // int mfindex = 0;
    // for(vector<Delta*>::iterator mf = (*f)->begin(); mf != (*f)->end(); mf++)
    // {
      // cout << "MF " << mfindex++ << " obs: " << (*mf)->observed() << " exp: " << (*mf)->expected() << " lamda: " << (*mf)->lambda() << endl;
    // }
  // }

  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    vector<int> featureXindices = _Xindices[(*f)->xListIndex()]; 
    vector<int> featureYindices = _Yindices[(*f)->yListIndex()]; 

    double xSize = 1.0;
    double ySize = 1.0;

    for(vector<int>::iterator x = featureXindices.begin(); x != featureXindices.end(); x++)
    {
      xSize *= _X->getBinSize(*x);
    }
    for(vector<int>::iterator y = featureYindices.begin(); y != featureYindices.end(); y++)
    {
      ySize *= _Y->getBinSize(*y);
    }

    (*f)->setRemainingAlphabetSize(xSize * ySize - 1.0);
    (*f)->setYAlphabetSize(ySize);
  }

}

int Model::nrOfFeatures()
{
  return features.size();
}

void Model::generateExpected()
{
  double sum = 0.0;
  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    for(vector<Delta*>::iterator mf = (*f)->begin(); mf != (*f)->end(); mf++)
    {
      sum += (*mf)->lambda();
    }
  }

  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    for(vector<Delta*>::iterator mf = (*f)->begin(); mf != (*f)->end(); mf++)
    {
      double e       = (*mf)->lambda();
      double zaehler = exp(e);
      int    xBar    = (*mf)->getUniqueXIndex();
      double count   = (*f)->getUniqueXCount(xBar);
      (*mf)->setExpected( count * zaehler / sum ); // zaehler / sum
    }
  }
}

ULContainer* Model::X()
{
  return _X;
}

ULContainer* Model::Y()
{
  return _Y;
}

ULContainer* Model::uniqueX(int i)
{
  return _uniqueX[i];
}

ULContainer* Model::uniqueY(int i)
{
  return _uniqueY[i];
}

Feature* Model::feature(int i)
{
  return features[i];
}
