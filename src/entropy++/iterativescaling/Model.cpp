#include <entropy++/iterativescaling/Model.h>

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
  // delete _uniqueX;
  // delete _uniqueY;
}



void Model::setRelations(vector<vector<int> > Xindices,
                         vector<vector<int> > Yindices,
                         vector<Relation>     relations)
{
  _Xindices  = Xindices;
  _Yindices  = Yindices;
  _relations = relations;
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
  for(int d = 0; d < _X->rows(); d++)
  {
    for(int r = 0; r < (int)_relations.size(); r++)
    {
      int i = _relations[r].first;
      int k = _relations[r].second;

      int j = _uniqueX[i]->find(_X, d);
      int l = _uniqueY[k]->find(_Y, d);

      bool found = false;
      for(vector<MFeature*>::iterator f = _features.begin();
          f != _features.end();
          f++)
      {
        if((*f)->match(i,j,k,l))
        {
          (*f)->incObserved();
          found = true;
          break;
        }
      }
      if(found == false)
      {
        MFeature *f = new MFeature(i,j,k,l);
        f->incObserved();
        _features.push_back(f);
      }
    }
  }
}

int Model::nrOfFeatures()
{
  return _features.size();
}
