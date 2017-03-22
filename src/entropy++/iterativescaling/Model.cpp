#include <entropy++/iterativescaling/Model.h>
// #include <omp.h>

// #include <glog/logging.h>

using namespace entropy;
using namespace entropy::iterativescaling;

# define EPSILON 0.00000001

Model::Model()
{
  _yAlphabetSize = 0;
  _conditional   = NULL;
  _marginal      = NULL;
  _joint         = NULL;
  _x_alphabet    = NULL;
  _y_alphabet    = NULL;

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
  delete Xalphabet;
  delete Yalphabet;

  delete _conditional;
  delete _marginal;
}

void Model::setFeatures(vector<vector<int> > Xindices,
    vector<vector<int> > Yindices, vector<Feature*> f)
{
  _Xindices = Xindices;
  _Yindices = Yindices;
  features = f;
}

void Model::createUniqueContainer()
{
 Xalphabet = Xdata->unique();
 Yalphabet = Ydata->unique();
  if (_x_alphabet != NULL) delete _x_alphabet;
  if (_y_alphabet != NULL) delete _y_alphabet;

  _x_indices.clear();
  _y_indices.clear();

  // get all columns for X
  for (vector<vector<int> >::iterator v = _Xindices.begin();
      v != _Xindices.end(); v++)
  {
    for (vector<int>::iterator i = (*v).begin(); i != (*v).end(); i++)
    {
      if (find(_x_indices.begin(), _x_indices.end(), *i) == _x_indices.end())
      {
        _x_indices.push_back(*i);
      }
    }
  }

  // get all columns for Y
  for (vector<vector<int> >::iterator v = _Yindices.begin();
      v != _Yindices.end(); v++)
  {
    for (vector<int>::iterator i = (*v).begin(); i != (*v).end(); i++)
    {
      if (find(_y_indices.begin(), _y_indices.end(), *i) == _y_indices.end())
      {
        _y_indices.push_back(*i);
      }
    }
  }
}

void Model::countObservedFeatures()
{

  for (int i = 0; i < Xdata->rows(); i++)
  {
    vector<unsigned long> xrow = Xdata->row(i);
    vector<unsigned long> yrow = Ydata->row(i);

    for (vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
    {
      bool found = false;
      for (vector<Delta*>::iterator d = (*f)->begin(); d != (*f)->end(); d++)
      {
        if ((*d)->matchXY(xrow, yrow))
        {
          (*d)->incObserved();
          found = true;
        }
      }
      if (found == false)
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

  deltaMatcher = new DeltaMatcher(Xdata->rows());
  for(int x = 0; x < Xdata->rows(); x++)
  {
    vector<unsigned long> x_row = Xdata->row(x);
    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchXY(x_row, y_row))
        {
          deltaMatcher->add(x, *d);
        }
      }
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

ULContainer* Model::XAlphabet()
{
  return Xalphabet;
}

ULContainer* Model::YAlphabet()
{
  return Yalphabet;
}

Feature* Model::feature(int i)
{
  return features[i];
}

void Model::calculateProbabilities()
{
  if (_conditional != NULL) delete _conditional;
  if (_marginal    != NULL) delete _marginal;
  if (_joint       != NULL) delete _joint;

  ULContainer *x_alphabet_full = Xdata->columns(_x_indices);
  _x_alphabet = x_alphabet_full->unique();

  ULContainer *y_alphabet_full = Ydata->columns(_y_indices);
  _y_alphabet = y_alphabet_full->unique();

  Matrix M(Xalphabet->rows(), Yalphabet->rows());

  for (int x = 0; x < Xalphabet->rows(); x++)
  {
    vector<unsigned long> x_row = Xalphabet->row(x);
    for (int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      for (vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if ((*d)->matchXY(x_row, y_row))
        {
          M(x,y) += (*d)->lambda();
        }
      }
    }
  }

  for (int x = 0; x < M.rows(); x++)
  {
    for (int y = 0; y < M.cols(); y++)
    {
      M(x, y) = exp(M(x,y));
    }
  }

  _conditional = new Matrix(_x_alphabet->rows(), _y_alphabet->rows());
  _joint       = new Matrix(_x_alphabet->rows(), _y_alphabet->rows());
  _marginal    = new Matrix(_x_alphabet->rows(), 1);

  for (int x = 0; x < Xalphabet->rows(); x++)
  {
    int x_row_index = __convertAlphabetToMatrixX(x);
    for (int y = 0; y < Yalphabet->rows(); y++)
    {
      int y_row_index = __convertAlphabetToMatrixY(y);
      (*_conditional)(x_row_index, y_row_index) += M(x, y);
    }
  }

  double sum = 0.0;
  for (int x = 0; x < _marginal->rows(); x++)
  {
    (*_marginal)(x, 0) = _conditional->rowSum(x);
    sum += (*_marginal)(x, 0);
    for (int y = 0; y < _conditional->cols(); y++)
    {
      (*_conditional)(x,y) = (*_conditional)(x,y) / (*_marginal)(x, 0);
    }
  }
  double b = x_alphabet_full->rows();
  for (int x = 0; x < _marginal->rows(); x++)
  {
    double a = x_alphabet_full->findlist(_x_alphabet, x).size();
    (*_marginal)(x, 0) = a / b;
  }

  delete x_alphabet_full;
  delete y_alphabet_full;
}

double Model::p_y_c_x(int yUniqueIndex, int xUniqueIndex)
{
  return (*_conditional)(xUniqueIndex, yUniqueIndex);
}

double Model::p_x(int xUniqueIndex)
{
  return (*_marginal)(xUniqueIndex, 0);
}

double Model::p_y_c_x_d(int yAlphabetIndex, int xAlphabetIndex)
{
  int x = __convertAlphabetToMatrixX(xAlphabetIndex);
  int y = __convertAlphabetToMatrixY(yAlphabetIndex);
  return (*_conditional)(x, y);
}

double Model::p_x_d(int xAlphabetIndex)
{
  int x = __convertAlphabetToMatrixX(xAlphabetIndex);
  return (*_marginal)(x, 0);
}

int Model::getNrOfUniqueX()
{
  return Xalphabet->rows();
}

int Model::getNrOfUniqueY()
{
  return Yalphabet->rows();
}

vector<int> Model::getAllColumnsForX()
{
  return _x_indices;
}

vector<int> Model::getAllColumnsForY()
{
  return _y_indices;
}

int Model::__convertAlphabetToMatrixX(int x)
{
  vector<unsigned long> x_row = Xalphabet->row(x);
  vector<unsigned long> x_values;
  for (vector<int>::iterator i = _x_indices.begin(); i != _x_indices.end(); i++)
  {
    x_values.push_back(x_row[*i]);
  }
  return _x_alphabet->find(x_values);
}

int Model::__convertAlphabetToMatrixY(int y)
{
  vector<unsigned long> y_row = Yalphabet->row(y);
  vector<unsigned long> y_values;
  for (vector<int>::iterator i = _y_indices.begin(); i != _y_indices.end(); i++)
  {
    y_values.push_back(y_row[*i]);
  }
  return _y_alphabet->find(y_values);
}

ULContainer* Model::__uniqueRows(vector<int> columns, ULContainer* data)
{
  vector<vector<unsigned long> > new_container;
  for(int r = 0; r < data->rows(); r++)
  {
    vector<unsigned long> row = data->row(r);
    bool found = false;
    for(vector<vector<unsigned long> >::iterator w = new_container.begin(); w != new_container.end(); w++)
    {
      bool rowEqual = true;
      for(int i = 0; i < (int)columns.size(); i++)
      {
        rowEqual &= ((*w)[columns[i]] == row[columns[i]]);
        if(rowEqual == false) break;
      }
      found |= rowEqual;
      if(found == true) break;
    }
    if(found == false)
    {
      new_container.push_back(row);
    }
  }
  ULContainer *nc = new ULContainer((int)new_container.size(), data->columns());

  for( vector<vector<unsigned long> >::iterator r = new_container.begin(); r != new_container.end(); r++)
  {
    for(vector<unsigned long>::iterator i = (*r).begin(); i != (*r).end(); i++)
    {
      *nc << *i;
    }
  }

  return nc;
}
