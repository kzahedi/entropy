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
  _s.resize(Yalphabet->rows());

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

    // cout << "X: ";
    // for(int i = 0; i < xrow.size() - 1; i++)
    // {
      // cout << xrow[i] << ", ";
    // }
    // cout << xrow[xrow.size()-1];
    // cout << "  --  Y: ";
    // for(int i = 0; i < yrow.size() - 1; i++)
    // {
      // cout << yrow[i] << ", ";
    // }
    // cout << yrow[yrow.size()-1];
    // cout << endl;


    for(vector<Feature*>::iterator f = features.begin(); f != features.end(); f++)
    {

      bool found = false;
      for(vector<Delta*>::iterator d = (*f)->begin(); d != (*f)->end(); d++)
      {
        if((*d)->matchXY(xrow,yrow))
        {
          (*d)->incObserved();
          // cout << "found: " << **d << endl;
          found = true;
        }
      }
      if(found == false)
      {
        vector<int> xIndices = _Xindices[(*f)->xListIndex()];
        vector<int> yIndices = _Yindices[(*f)->yListIndex()];
        Delta *d = new Delta(xrow, xIndices, yrow, yIndices);
        d->incObserved();
        // cout << "adding: " << *d << endl;
        deltas.push_back(d);
        (*f)->push_back(d);
      }
    }
  }

  // for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  // {
    // (*d)->setObserved((*d)->observed() / (double)Xdata->rows());
  // }
}

int Model::nrOfFeatures()
{
  return features.size();
}

void Model::generateExpected()
{

  for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
  {
    (*d)->setExpected(0.0);
  }

  for(int j = 0; j < Xdata->rows(); j++)
  {
    vector<unsigned long> x_row = Xdata->row(j);

    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      _s[y] = 0.0;
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchXY(x_row, y_row))
        {
          _s[y] += (*d)->lambda();
        }
      }
    } // for each output y

    // double z = _yAlphabetSize - (int)_s.size();
    double z = 0.0;
    for(vector<double>::iterator i = _s.begin(); i != _s.end(); i++) z += exp(*i);

    for(int y = 0; y < Yalphabet->rows(); y++)
    {
      vector<unsigned long> y_row = Yalphabet->row(y);
      for(vector<Delta*>::iterator d = deltas.begin(); d != deltas.end(); d++)
      {
        if((*d)->matchXY(x_row, y_row))
        {
          (*d)->setExpected((*d)->expected() + exp(_s[y]) / z);
        }
      }
    }
  } // j
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

