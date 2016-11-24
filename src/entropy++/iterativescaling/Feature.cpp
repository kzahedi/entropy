#include "Feature.h"

using namespace entropy::iterativescaling;

Feature::Feature()
{
  _lambda = new SparseMatrix();
  _sizeX  = 0;
  _sizeY  = 0;
}

//(Eingabealphabet, Anzahl der X und Y Werte, Anzahl der Testwerte, ein Startwert fuer alle  lambda)
Feature::Feature(int sizeX, int sizeY, double valuelambda) // systXsize, int systYsize, double valuelambda)
{
  _sizeX  = sizeX;
  _sizeY  = sizeY;
  _lambda = new SparseMatrix(valuelambda);
}

//alle lambda explizit ueber die Matrix setzen
Feature::Feature(int sizeX, int sizeY, SparseMatrix &lambda) // int systXsize, int systYsize, SparseMatrix &lambda)
{
  _sizeX  = sizeX; // pow(xAlphabet.rows(),systXsize);
  _sizeY  = sizeY; // pow(yAlphabet.rows(),systYsize);
  _lambda = &lambda;
}

Feature::~Feature()
{
  delete _lambda;
}

int Feature::getLambdaSize()
{
  int j = _lambda->size();
  return j;
}

double Feature::getLambda(int i, int j)
{
  assert(i < _sizeX && j < _sizeY);
  return (*_lambda)(i,j);
}

void Feature::setLambda(int i, int j, double newvalue)
{
  assert(i < _sizeX && j < _sizeY);
  (*_lambda)(i,j) = newvalue;
}

Feature& Feature::operator=(const Feature& c)
{
  _sizeX = c._sizeX;
  _sizeY = c._sizeY;

  _lambda     = new SparseMatrix();
  _lambda     = c._lambda;

  return *this;
}

