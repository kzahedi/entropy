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

  Xdata = X;
  Ydata = Y;

}

Model::~Model()
{
  delete Xdata;
  delete Ydata;

  delete Xalphabet;
  delete Yalphabet;
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
  Xalphabet = Xdata->unique();
  Yalphabet = Ydata->unique();
}

void Model::countObservedFeatures()
{
  for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
  {
    int i = (*f)->xListIndex();
    int k = (*f)->yListIndex();

    vector<int> xIndices = _Xindices[i];
    vector<int> yIndices = _Yindices[k];

    for(vector<int>::iterator xx = xIndices.begin(); xx != xIndices.end(); xx++)
    {
      cout << *xx << " ";
    }
    cout << endl;

    for(int d = 0; d < Xdata->rows(); d++)
    {
      vector<unsigned long> xrow = Xdata->row(d);
      vector<unsigned long> yrow = Ydata->row(d);

      bool found = false;
      for(vector<Delta*>::iterator d = (*f)->begin(); d != (*f)->end(); d++)
      {
        if((*d)->matchXY(xrow,yrow))
        {
          (*d)->incObserved();
          found = true;
          break;
        }
      }
      if(found == false)
      {
        Delta *d = new Delta(xrow, xIndices, yrow, yIndices);
        d->incObserved();
        (*f)->push_back(d);
        deltas.push_back(d);
      }
    }
  }
}

int Model::nrOfFeatures()
{
  return features.size();
}

void Model::generateExpected()
{
  for(int j = 0; j < Xdata->rows(); j++)
  {
    double z = 0.0;

    vector<unsigned long> x_row = Xdata->row(j);

    for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
    {
      if((*d)->matchX(x_row))
      {
        z += (*d)->lambda();
      }
    }

    z = exp(z);
    for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
    {
      if((*d)->matchX(x_row))
      {
        double e = (*d)->expected();
        // cout << "1. e: " << e << endl;
        e += (*d)->lambda() / z;
        // cout << "2. e: " << e << endl;
        (*d)->setExpected(e);
      }
    }
  }
}


ULContainer* Model::X()
{
  return Xdata;
}

ULContainer* Model::Y()
{
  return Ydata;
}

Feature* Model::feature(int i)
{
  return features[i];
}

void Model::calculateProbabilities()
{
  // _conditionals = new Matrix(XalphabetFromData->rows(), YalphabetFromData->rows()); // TODO: clean up in destructor
  // _marginals    = new Matrix(XalphabetFromData->rows(), 1);

  // double nenner = 0.0;

  // for(int x = 0; x < XalphabetFromData->rows(); x++)
  // {
    // for(int y = 0; y < YalphabetFromData->rows(); y++)
    // {
      // double sum = 0.0;
      // for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
      // {
        // int xListIndex = (*f)->xListIndex();
        // int yListIndex = (*f)->yListIndex();

        // vector<int> fx = Xalphabet[xListIndex]->findlist(XalphabetFromDataPerFeature[xListIndex], x);
        // vector<int> fy = Yalphabet[yListIndex]->findlist(YalphabetFromDataPerFeature[yListIndex], y);

        // for(vector<int>::iterator i = fx.begin(); i != fx.end(); i++)
        // {
          // for(vector<int>::iterator j = fy.begin(); j != fy.end(); j++)
          // {
            // for(vector<Delta*>::iterator d = (*f)->begin(); d != (*f)->end(); d++)
            // {
              // if((*d)->match(*i, *j))
              // {
                // sum += (*d)->lambda();
              // }
            // }
          // }
        // }
      // }
      // if(fabs(sum) <= EPSILON)
      // {
        // (*_conditionals)(x,y) = 0.0; // TODO check
      // }
      // else
      // {
        // (*_conditionals)(x,y) = exp(sum);
      // }
    // }
  // }

  // double sum = 0.0;
  // for(int x = 0; x < _conditionals->rows(); x++)
  // {
    // (*_marginals)(x,0) = _conditionals->rowSum(x);
    // sum += (*_marginals)(x,0);
    // for(int y = 0; y < _conditionals->cols(); y++)
    // {
      // (*_conditionals)(x,y) = (*_conditionals)(x,y) / (*_marginals)(x,0);
    // }
  // }

  // for(int x = 0; x < _marginals->rows(); x++)
  // {
    // (*_marginals)(x,0) = (*_marginals)(x,0) / sum;
  // }
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
  return Xalphabet->rows();
}

int Model::getNrOfUniqueY()
{
  return Yalphabet->rows();
}

