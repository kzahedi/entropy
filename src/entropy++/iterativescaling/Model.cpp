#include <entropy++/iterativescaling/Model.h>

using namespace entropy::iterativescaling;

// obs(x)
void Feature::setUniqueXCount(int index, int count)
{
  if(_uniqueXCount.size() < index + 1)
  {
    _uniqueXCount.resize(index+1);
  }
  _uniqueXCount[index] = count;
}

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
  // delete _uniqueX;
  // delete _uniqueY;
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
      for(vector<MFeature*>::iterator mf = (*f)->begin(); mf != (*f)->end(); mf++)
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
        MFeature *mf = new MFeature(j,l);
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

  int index = 0;
  cout << "In countObservedFeatures: "<< endl;
  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    cout << "Feature " << index++ << endl;
    int i = (*f)->xListIndex(); // list of columns that define X
    for(int x = 0; x < _uniqueX[i]->rows(); x++)
    {
      cout << "  " << x << " x count: " << (*f)->getUniqueXCount(x) << endl;
    }
    int mfindex = 0;
    for(vector<MFeature*>::iterator mf = (*f)->begin(); mf != (*f)->end(); mf++)
    {
      cout << "MF " << mfindex++ << " obs: " << (*mf)->observed() << " exp: " << (*mf)->expected() << " lamda: " << (*mf)->lambda() << endl;
    }
  }

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

    (*f)->setRemainingAlphabetSize(xSize * ySize - _uniqueX[(*f)->xListIndex()]->rows());
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
    for(vector<MFeature*>::iterator mf = (*f)->begin(); mf != (*f)->end(); mf++)
    {
      double e       = (*mf)->lambda() - (*f)->getRemainingAlphabetSize();
      double zaehler = exp(e);
      (*mf)->setExpected(zaehler);
      sum += zaehler;
    }
  }

  cout << "Sum: " << sum << endl;

  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    for(vector<MFeature*>::iterator mf = (*f)->begin();
        mf != (*f)->end(); mf++)
    {
      int    xBar  = (*mf)->getUniqueXIndex();
      double count = (*f)->getUniqueXCount(xBar);
      (*mf)->setExpected( count * (*mf)->expected() / sum ); // zaehler / sum
    }
  }

  int index = 0;
  cout << "In generateExpected: "<< endl;
  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    cout << "Feature " << index++ << endl;
    int mfindex = 0;
    for(vector<MFeature*>::iterator mf = (*f)->begin(); mf != (*f)->end(); mf++)
    {
      cout << "MF " << mfindex++ << " obs: " << (*mf)->observed() << " exp: " << (*mf)->expected() << " lamda: " << (*mf)->lambda() << endl;
    }
  }

}
