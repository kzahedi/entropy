#include "IT.h"

// IT::IT(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,double lambdavalue, bool gis)
IT::IT(DContainer &eX, DContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param, bool gis)
{
  assert(eX.rows() == eY.rows());
  _param       = param;
  _gis         = gis;
  _valX        = &eX;
  _valY        = &eY;
  _X           = &aX;
  _Y           = &aY;
  _sizeSystX   = systX.size();
  _sizeX       = _X->rows();
  _sizeY       = _Y->rows();
  _sizeColValY = _valY->columns();
  _sizeColValX = _valX->columns();
  _sizeRowValX = _valX->rows();
  _sizeRowValY = _valY->rows();
  _systX       = systX;
  _systY       = systY;

  if(_gis) // gis and csgis require different feature matrices
  {
    _FM = new FeatureMatrix(*_valX,*_valY,*_X,*_Y, systX, systY, param.lambdavalue);
  }
  else
  {
    _IM = new InstanceMatrix(*_valX,*_valY,*_X,*_Y, systX, systY, param.lambdavalue);
  }
  _observed=__getobs();
}
// Konstruktor fuer Alphabete mit unsigned long
IT::IT(ULContainer &eX, ULContainer &eY, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY, IsParameter param, bool gis)
{
  assert(eX.rows() == eY.rows());
  _param       = param;
  _gis         = gis;
  _valX        = NULL;
  _valY        = NULL;
  ULContainer *valX = &eX;
  ULContainer *valY = &eY;
  _X           = &aX;
  _Y           = &aY;
  _sizeSystX   = systX.size();
  _sizeX       = _X->rows();
  _sizeY       = _Y->rows();
  _sizeColValY = valY->columns();
  _sizeColValX = valX->columns();
  _sizeRowValX = valX->rows();
  _sizeRowValY = valY->rows();
  _systX       = systX;
  _systY       = systY;
  if(_gis) // gis and csgis require different feature matrices
  {
    _FM = new FeatureMatrix(*valX,*valY,*_X,*_Y, systX, systY, param.lambdavalue);
  }
  else
  {
    _IM = new InstanceMatrix(*valX,*valY,*_X,*_Y, systX, systY, param.lambdavalue);
  }
  _observed=__getobs();
}

IT::IT(int ColValY, DContainer &eX, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY)
{
  _systX       = systX;
  _systY       = systY;
  _gis         = true;
  _X           = &aX;
  _Y           = &aY;
  _sizeSystX   = systX.size();
  _sizeX       = _X->rows();
  _sizeY       = _Y->rows();
  _valX        = &eX;
  _sizeColValY = ColValY;
  _sizeColValX = (*_valX).columns();
  _sizeRowValX = (*_valX).rows();
  _valY        = new DContainer(_sizeRowValX,ColValY);
  _sizeRowValY = _sizeRowValX;
  _FM          = new FeatureMatrix(*_valX,*_valY,*_X,*_Y,systX,systY,1);
  _IM          = NULL;
  _observed    = NULL;
}
// p(y | x)
// double IT::prop(int rowX, vector<vector<double> > Y, int rowY)
double IT::prop(int rowX, int rowY)
{
  double featexp   = 0;
  double featnorm = 0;
  double norm     = 0;
  double exponent = 0;
  for(int feat=0; feat< _systX.size(); feat++)
  {
      if(_gis)
      {
        featexp += (*_FM).getFeatureArrayvalueAlphY(feat,rowX, rowY);
      }
      else
      {
        featexp += (*_IM).getFeatureArrayvalueAlphY(feat,rowX, rowY);
      }
  }
  exponent = exp(featexp);
  for(int yi=0;yi<pow(_Y->rows(),_sizeColValY);yi++)
  {
    for(int feat=0; feat< _systX.size(); feat++)
    {
        if(_gis)
        {
          featnorm += (*_FM).getFeatureArrayvalueAlphY(feat,rowX,yi);
        }
        else
        {
          featnorm += (*_IM).getFeatureArrayvalueAlphY(feat,rowX,yi);
        }
    }
    norm     += exp(featnorm);
    featnorm  = 0;
  }
  return exponent/norm;
}
// p(y | x)
// double IT::prop(int indexX, vector<vector<double> > Y, int indexY)
double IT::propAlphX(int indexX, int rowY)
{
  double featexp   = 0;
  double featnorm = 0;
  double norm     = 0;
  double exponent = 0;
  for(int feat=0; feat< _systX.size(); feat++)
  {
      if(_gis)
      {
        featexp += (*_FM).getFeatureArrayvalueAlphYAlphX(feat,indexX, rowY);
      }
      else
      {
        featexp += (*_IM).getFeatureArrayvalueAlphYAlphX(feat,indexX, rowY);
      }
  }
  exponent = exp(featexp);
  for(int yi=0;yi<pow(_Y->rows(),_sizeColValY);yi++)
  {
    for(int feat=0; feat< _systX.size(); feat++)
    {
        if(_gis)
        {
          featnorm += (*_FM).getFeatureArrayvalueAlphYAlphX(feat,indexX,yi);
        }
        else
        {
          featnorm += (*_IM).getFeatureArrayvalueAlphYAlphX(feat,indexX,yi);
        }
    }
    norm     += exp(featnorm);
    featnorm  = 0;
  }
  return exponent/norm;
}
// P(x)
double IT::propm(int rowX){
  double z=0;
  double exponent=0;
  for(int y=0;y<pow(_Y->rows(),_sizeColValY);y++)
  {
    for(int feat=0;feat<  _systX.size();feat++)
    {
        if(_gis)
        {
          exponent+=(*_FM).getFeatureArrayvalueAlphYAlphX(feat,rowX,y);
        }
        else
        {
          exponent+=(*_IM).getFeatureArrayvalueAlphYAlphX(feat,rowX,y);
        }

    }
    z+=exp(exponent);
    exponent=0;
  }
  double n;
  double nexp=0;
  for(int x=0;x<pow(_X->rows(),_sizeColValX);x++)
  {
    for(int y=0;y<pow(_Y->rows(),_sizeColValY);y++)
    {
        for(int feat=0;feat< _systX.size();feat++)
        {
          if(_gis)
          {
            exponent+=(*_FM).getFeatureArrayvalueAlphYAlphX(feat,x,y);
          }
          else
          {
            exponent+=(*_IM).getFeatureArrayvalueAlphYAlphX(feat,x,y);
          }

      }
      nexp +=exp(exponent);
      exponent=0;
    }
    n+=nexp;
    nexp=0;
  }
  return z/n;
}

double IT::getFeatureArraylambda(int Feati, int ilambdaX, int ilambdaY)
{
  assert(Feati<_systX.size()); //kontrolle der Parameter ilambdaX und ilambdaY findet in Feature statt
  double lambda = 0.0;
  if(_gis)
  {
    lambda = _FM->getFeatureArraylambda(Feati, ilambdaX, ilambdaY);
  }
  else
  {
    lambda = _IM->getFeatureArraylambda(Feati, ilambdaX, ilambdaY);
  }
  return lambda;
}

// get observed
double*** IT::__getobs()
{
  //erzeugen von observed
  _observed = new double**[_sizeSystX];
  for(int i=0; i<_systX.size(); i++)
  {
      _observed[i]=new double*[(int) pow(_X->rows(),_systX[i].size())];
      for(int k=0; k< pow(_X->rows(),_systX[i].size()); k++)
      {
        _observed[i][k]= new double[(int) pow(_Y->rows(),_systY[i].size())];
        for(int l=0; l< pow(_Y->rows(),_systY[i].size());l++)
        {
          _observed[i][k][l] = 0.0;
        }
      }
   }
  // fuellen von observed
  for(int i=0;i<_sizeRowValX;i++ )
  {
    for(int feat=0; feat< _sizeSystX;feat++)
    {
        for(int delti = 0; delti<pow(_X->rows(),_systX[feat].size()); delti++)
        {
          for(int deltj=0; deltj<pow(_Y->rows(),_systY[feat].size()); deltj++)
          {
            if(_gis)
            {
              if(_FM->getFeatureArraydelta(feat,delti,deltj,i,i)==1)
              {
                _observed[feat][delti][deltj]++;
              }
            }
            else
            {
              if(_IM->getFeatureArraydelta(feat,delti,deltj,i,i)==1)
              {
                _observed[feat][delti][deltj]++;
              }
            }
          }
        }
      }
  }
  return _observed;
}
IT::~IT(){
	if(_observed!=NULL)
	  {
	    for(int i=0;i<_sizeSystX;i++)
	    {
	      for(int k=0;k<pow(_X->rows(),_systX[i].size());k++)
	      {
	        delete [] _observed[i][k];
	      }
	      delete [] _observed[i];
	    }
	    delete [] _observed;
	  }
	  for(int i=0; i< _systX.size();i++)
	  {
	    _systX[i].clear();
	  }
	  _systX.clear();
	  for(int j=0; j<_systY.size();j++)
	  {
		_systY[j].clear();
	  }
	  _systY.clear();
}
//veraendern der lambdas im Algorithmus funktioniert mit einer Methode in ITMatrix, hier nur zum erzeugen von Testwerten
void IT::setFeatureArraylambda(int feati, int ilambdaX, int ilambdaY,double valuelambda){
  assert(feati<_systX.size());
  assert(_gis == true);
  _FM->setFeatureArraylambda(feati,ilambdaX,ilambdaY,valuelambda);
}

// Berechnung der Alphabetreihe aus dem Index
vector<double> IT::index(int index,bool x,int sizeCol)
{
  int sizeAlph;
  vector<double> zeile;
  double z;
  if(x){
	  assert(index<pow(_X->rows(),sizeCol));
	  sizeAlph=_sizeX;
  }
  else{
	  assert(index<pow(_Y->rows(),sizeCol));
	  sizeAlph=_sizeY;
  }
  for(int i = sizeCol; i>0 ; i--)
  {
	  z = index % (int) (pow(sizeAlph,i));
	  z = z/ (pow(sizeAlph,i-1));
	  int j= (int) z;
	  if(x){
		  zeile.push_back((*_X)(j,0));
	  }
	  else{
		  zeile.push_back((*_Y)(j,0));
	  }
  }
  return zeile;
}
