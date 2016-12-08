#include <entropy++/iterativescaling/Model.h>

using namespace entropy::iterativescaling;

void Feature::setUniqueXCount(int index, int count)
{
  if(_uniqueXCount.size() < index)
  {
    _uniqueXCount.resize(index);
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
                        vector<Feature*>     features)
{
  _Xindices = Xindices;
  _Yindices = Yindices;
  _features = features;
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
  for(vector<Feature*>::iterator f = _features.begin(); f != _features.end(); f++)
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
  for(vector<Feature*>::iterator f = _features.begin(); f != _features.end(); f++)
  {
    int i = (*f)->xListIndex(); // list of columns that define X
    ULContainer* xUniqueForThisFeature = _uniqueX[i];
    ULContainer* xMatchColumns         = _X->columns(_Xindices[i]);
    int count = xMatchColumns->findlist(xUniqueForThisFeature,0).size();
    (*f)->setUniqueXCount(i, count);
  }

}

int Model::nrOfFeatures()
{
  return features.size();
}

void Model::generateExpected()
{
  double sum = 0.0;
  for(vector<Feature*>::iterator f = _features.begin(); f != _features.end(); f++)
  {
    double sum_of_lambdas = 0.0;
    for(vector<MFeature*>::iterator mf = (*f)->begin(); mf != (*f)->end(); mf++)
    {
      int xBar     = (*mf)->getUniqueXIndex();
      double count = (*f)->getUniqueXCount(xBar);
      sum_of_lambdas += count * (*mf)->lambda();
      double zaehler = exp(sum_of_lambdas);
      (*mf)->setExpected(zaehler);
      sum += zaehler;
    }

    for(vector<MFeature*>::iterator f = features.begin();
        f != features.end(); f++)
    {
      (*f)->setExpected( (*f)->expected() / sum ); // zaehler / sum
    }
  }
}
