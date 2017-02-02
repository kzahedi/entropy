#include <entropy++/iterativescaling/Model.h>

#include <glog/logging.h>

using namespace entropy;
using namespace entropy::iterativescaling;

# define EPSILON 0.00000001

Model::Model()
{
  _yAlphabetSize = 0;
  _conditionals  = NULL;
  _marginals     = NULL;
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
  // delete Xdata;
  // delete Ydata;

  delete Xalphabet;
  delete Yalphabet;

  delete _conditionals;
  delete _marginals;
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

  _yAlphabetSize = 1.0;
  for(int c = 0; c < Yalphabet->columns(); c++)
  {
    _yAlphabetSize *= Yalphabet->getBinSize(c);
  }
}

void Model::countObservedFeatures()
{
  for(int d = 0; d < Xdata->rows(); d++)
  {
    vector<unsigned long> xrow = Xdata->row(d);
    vector<unsigned long> yrow = Ydata->row(d);

    for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
    {
      bool found = false;
      for(vector<Delta*>::iterator d = (*f)->begin(); d != (*f)->end(); d++)
      {
        if((*d)->matchXY(xrow,yrow))
        {
          (*d)->incObserved();
          found = true;
        }
      }
      if(found == false)
      {
        vector<int> xIndices = _Xindices[(*f)->xListIndex()];
        vector<int> yIndices = _Yindices[(*f)->yListIndex()];
        Delta *d = new Delta(xrow, xIndices, yrow, yIndices);
        d->incObserved();
        deltas.push_back(d);
        (*f)->push_back(d);
      }
    }
  }
  if(VLOG_IS_ON(100))
  {
    VLOG(100) << "Counting done";
    for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
    {
      VLOG(100) << **d;
    }
  }
}

int Model::nrOfFeatures()
{
  return features.size();
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
  if(_conditionals != NULL) delete _conditionals;
  if(_marginals    != NULL) delete _marginals;

  // get all columns for X
  vector<int> x_indices;
  for(vector<vector<int> >::iterator v = _Xindices.begin(); v != _Xindices.end(); v++)
  {
    for(vector<int>::iterator i = (*v).begin(); i != (*v).end(); i++)
    {
      if(find(x_indices.begin(), x_indices.end(), *i) == x_indices.end())
      {
        x_indices.push_back(*i);
      }
    }
  }

  // get all columns for Y
  vector<int> y_indices;
  for(vector<vector<int> >::iterator v = _Yindices.begin(); v != _Yindices.end(); v++)
  {
    for(vector<int>::iterator i = (*v).begin(); i != (*v).end(); i++)
    {
      if(find(y_indices.begin(), y_indices.end(), *i) == y_indices.end())
      {
        y_indices.push_back(*i);
      }
    }
  }

  ULContainer *x_alphabet_full = Xalphabet->columns(x_indices);
  ULContainer *y_alphabet_full = Yalphabet->columns(y_indices);
  ULContainer *x_alphabet = x_alphabet_full->unique();
  ULContainer *y_alphabet = y_alphabet_full->unique();

  Matrix M(Xalphabet->rows(), Yalphabet->rows());

  for(int x = 0; x < Xalphabet->rows(); x++)
  {
    vector<unsigned long> x_row = Xalphabet->row(x);
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      double sum = 0.0;
      for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchXY(x_row, y_row))
        {
          sum += (*d)->lambda();
        }
      }
      if(fabs(sum) <= EPSILON)
      {
        M(x,y) = 0.0; // TODO check
      }
      else
      {
        M(x,y) = exp(sum);
      }
    }
  }

  _conditionals = new Matrix(x_alphabet->rows(), y_alphabet->rows());
  _marginals    = new Matrix(x_alphabet->rows(), 1);

  for(int x = 0; x < Xalphabet->rows(); x++)
  {
    vector<unsigned long> x_row = Xalphabet->row(x);
    vector<unsigned long> x_values;
    for(vector<int>::iterator i = x_indices.begin(); i != x_indices.end(); i++)
    {
      x_values.push_back(x_row[*i]);
    }
    int x_row_index = x_alphabet->find(x_values);
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      vector<unsigned long> y_values;
      for(vector<int>::iterator i = y_indices.begin(); i != y_indices.end(); i++)
      {
        y_values.push_back(y_row[*i]);
      }
      int y_row_index = y_alphabet->find(y_values);
      (*_conditionals)(x_row_index, y_row_index) += M(x,y);
    }
  }

  double sum = 0.0;
  for(int x = 0; x < _marginals->rows(); x++)
  {
    (*_marginals)(x,0) = _conditionals->rowSum(x);
    sum += (*_marginals)(x,0);
    for(int y = 0; y < _conditionals->cols(); y++)
    {
      (*_conditionals)(x,y) = (*_conditionals)(x,y) / (*_marginals)(x,0);
    }
  }

  double b = x_alphabet_full->rows();
  for(int x = 0; x < _marginals->rows(); x++)
  {
    double a = x_alphabet_full->findlist(x_alphabet, x).size();
    (*_marginals)(x,0) = a / b;
  }

  delete x_alphabet;
  delete y_alphabet;
  delete x_alphabet_full;
  delete y_alphabet_full;

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

