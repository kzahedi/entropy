#include <entropy++/iterativescaling/Model.h>

using namespace entropy;
using namespace entropy::iterativescaling;

# define EPSILON 0.00000001

Model::Model()
{
  _yAlphabetSize = 0;
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


  _uniqueXFromData = _X->unique();
  _uniqueYFromData = _Y->unique();

  _uniqueXFromDataPerFeature = new ULContainer*[_Xindices.size()];
  _uniqueYFromDataPerFeature = new ULContainer*[_Yindices.size()];

  _xFromDataPerFeature = new ULContainer*[_Xindices.size()];
  _yFromDataPerFeature = new ULContainer*[_Yindices.size()];

  for(int i = 0; i < _Xindices.size(); i++)
  {
    _uniqueXFromDataPerFeature[i] = _uniqueXFromData->columns(_Xindices[i]);
  }
  for(int i = 0; i < _Yindices.size(); i++)
  {
    _uniqueYFromDataPerFeature[i] = _uniqueYFromData->columns(_Yindices[i]);
  }

  for(int i = 0; i < _Xindices.size(); i++)
  {
    _xFromDataPerFeature[i] = _X->columns(_Xindices[i]);
  }
  for(int i = 0; i < _Yindices.size(); i++)
  {
    _yFromDataPerFeature[i] = _Y->columns(_Yindices[i]);
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

    (*f)->setXAlphabetSize(xSize);
    (*f)->setYAlphabetSize(ySize);
  }

  _yAlphabetSize = 1;
  for(int c = 0; c < _Y->columns(); c++)
  {
    _yAlphabetSize *= _Y->getBinSize(c);
  }

}

int Model::nrOfFeatures()
{
  return features.size();
}

void Model::generateExpected()
{
  for(int j = 0; j < _uniqueXFromData->rows(); j++)
  {
    double z = 0.0;
    double y_count = 0.0;
    vector<int> row_indices = _X->findlist(_uniqueXFromData, j); // TODO do this only once

    for(vector<int>::iterator i = row_indices.begin(); i != row_indices.end(); i++)
    {
      for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
      {
        int xListIndex = (*f)->xListIndex();
        int yListIndex = (*f)->yListIndex();

        vector<int> fx = _uniqueX[xListIndex]->findlist(_xFromDataPerFeature[xListIndex],  j);
        vector<int> fy = _uniqueY[yListIndex]->findlist(_yFromDataPerFeature[yListIndex], *i);

        for(vector<int>::iterator fi = fx.begin(); fi != fx.end(); fi++)
        {
          for(vector<int>::iterator fj = fy.begin(); fj != fy.end(); fj++)
          {
            for(vector<Delta*>::iterator d = (*f)->begin(); d != (*f)->end(); d++)
            {
              if((*d)->match(*fi, *fj))
              {
                z += (*d)->lambda();
                y_count++;
              }
            }
          }
        }
      }
    }

    z = exp(z) + (_yAlphabetSize - y_count);

    for(vector<int>::iterator i = row_indices.begin(); i != row_indices.end(); i++)
    {
      for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
      {
        int xListIndex = (*f)->xListIndex();
        int yListIndex = (*f)->yListIndex();

        vector<int> fx = _uniqueX[xListIndex]->findlist(_xFromDataPerFeature[xListIndex],  j);
        vector<int> fy = _uniqueY[yListIndex]->findlist(_yFromDataPerFeature[yListIndex], *i);

        for(vector<int>::iterator fi = fx.begin(); fi != fx.end(); fi++)
        {
          for(vector<int>::iterator fj = fy.begin(); fj != fy.end(); fj++)
          {
            for(vector<Delta*>::iterator d = (*f)->begin(); d != (*f)->end(); d++)
            {
              if((*d)->match(*fi, *fj))
              {
                int    xBar  = (*d)->getUniqueXIndex();
                double count = (*f)->getUniqueXCount(xBar);
                double e     = (*d)->expected();
                e += exp((*d)->lambda()) / z;
                (*d)->setExpected(e);
              }
            }
          }
        }
      }
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

void Model::calculateProbabilities()
{
  _conditionals = new Matrix(_uniqueXFromData->rows(), _uniqueYFromData->rows()); // TODO: clean up in destructor
  _marginals    = new Matrix(_uniqueXFromData->rows(), 1);

  double nenner = 0.0;

  for(int x = 0; x < _uniqueXFromData->rows(); x++)
  {
    for(int y = 0; y < _uniqueYFromData->rows(); y++)
    {
      double sum = 0.0;
      for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
      {
        int xListIndex = (*f)->xListIndex();
        int yListIndex = (*f)->yListIndex();

        vector<int> fx = _uniqueX[xListIndex]->findlist(_uniqueXFromDataPerFeature[xListIndex], x);
        vector<int> fy = _uniqueY[yListIndex]->findlist(_uniqueYFromDataPerFeature[yListIndex], y);

        for(vector<int>::iterator i = fx.begin(); i != fx.end(); i++)
        {
          for(vector<int>::iterator j = fy.begin(); j != fy.end(); j++)
          {
            for(vector<Delta*>::iterator d = (*f)->begin(); d != (*f)->end(); d++)
            {
              if((*d)->match(*i, *j))
              {
                sum += (*d)->lambda();
              }
            }
          }
        }
      }
      if(fabs(sum) <= EPSILON)
      {
        (*_conditionals)(x,y) = 0.0; // TODO check
      }
      else
      {
        (*_conditionals)(x,y) = exp(sum);
      }
    }
  }

  double sum = 0.0;
  for(int x = 0; x < _conditionals->rows(); x++)
  {
    (*_marginals)(x,0) = _conditionals->rowSum(x);
    sum += (*_marginals)(x,0);
    for(int y = 0; y < _conditionals->cols(); y++)
    {
      (*_conditionals)(x,y) = (*_conditionals)(x,y) / (*_marginals)(x,0);
    }
  }

  for(int x = 0; x < _marginals->rows(); x++)
  {
    (*_marginals)(x,0) = (*_marginals)(x,0) / sum;
  }
}

double Model::p_y_c_x(int yUniqueIndex, int xUniqueIndex)
{
  return (*_conditionals)(xUniqueIndex, yUniqueIndex);
}

double Model::p_x(int xUniqueIndex)
{
  return (*_marginals)(xUniqueIndex,0);
}

int Model::getNrOfUniqueX()
{
  return _uniqueXFromData->rows();
}

int Model::getNrOfUniqueY()
{
  return _uniqueYFromData->rows();
}

