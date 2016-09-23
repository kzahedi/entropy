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
   // _IM = new InstanceMatrix(*_valX,*_valY,*_X,*_Y, param.lambdavalue);
  }
  _observed=__getobs();
}
//umstellen
IT::IT(int ColValY, DContainer &eX, DContainer &aX, DContainer &aY,vector<vector<int> > systX, vector<vector<int> > systY)
{
  _systX       = systX;
  _systY       = systY;
  _gis         = true;
  _X           = &aX;
  _Y           = &aY;
  _sizeX       = _X->rows();
  _sizeY       = _Y->rows();
  _valX        = &eX;
  _sizeColValY = ColValY;
  _sizeColValX = (*_valX).columns();
  _sizeRowValX = (*_valX).rows();
  _valY        = new DContainer(_sizeRowValX,ColValY);
  _sizeRowValY = 0;
  _FM          = new FeatureMatrix(*_valX,*_valY,*_X,*_Y,systX,systY,1);
  _IM          = NULL;
  _observed    = NULL;
}
/*
// returns p(y_j = valY | x_i = valX)
double IT::prop(int Feati, int Featj, double ValX, double ValY)
{
  assert(Feati<_sizeColValX && Featj<_sizeColValY);
  double norm     = 0;
  double exponent = 0;
  if(_gis)
  {
    exponent = exp((*_FM).getFeatureArrayvalue(Feati,Featj,ValX,ValY));
  }
  else
  {
    exponent = exp((*_IM).getFeatureArrayvalue(Feati,Featj,ValX,ValY));
  }

  for(int yi = 0; yi < _sizeY; yi++)
  {
    if(_gis)
    {
      norm += exp((*_FM).getFeatureArrayvalue(Feati,Featj,ValX,(*_Y)(yi,0)));
    }
    else
    {
      norm += exp((*_IM).getFeatureArrayvalue(Feati,Featj,ValX,(*_Y)(yi,0)));
    }
  }
  return exponent/norm;
} */

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
     //   featexp += (*_IM).getFeatureArrayvalueAlphY(feat,rowX, rowY);
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
    //      featnorm += (*_IM).getFeatureArrayvalue(feat,rowX,yi);
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
     //     exponent+=(*_IM).getFeatureArrayvalueAlphYAlphX(feat,rowX,y);
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
     //       exponent+=(*_IM).getFeatureArrayvalue(feat,x,y);
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
  assert(Feati<_systX.size());
  double lambda = 0.0;
  if(_gis)
  {
    lambda = _FM->getFeatureArraylambda(Feati, ilambdaX, ilambdaY);
  }
  else
  {
 //   lambda = _IM->getFeatureArraylambda(Feati, ilambdaX, ilambdaY);
  }
  return lambda;
}

// get observed
double*** IT::__getobs()
{
  _observed = new double**[_systX.size()];
  for(int i=0; i<_systX.size(); i++)
  {
      _observed[i]=new double*[_sizeX];
      for(int k=0; k< _sizeX; k++)
      {
        _observed[i][k]= new double[_sizeY];
        for(int l=0; l< _sizeY;l++)
        {
          _observed[i][k][l] = 0.0;
        }
      }
   }
  vector<double> x;
  vector<double> y;
  //vector observed
  for(int i=0;i<_sizeRowValX;i++ )
  {
    for(int feat=0; feat< _systX.size();feat++)
    {
        for(int delti = 0; delti < pow(_X->rows(),_systX[feat].size()); delti++)
        {
          for(int deltj=0; deltj<pow(_Y->rows(),_systY[feat].size()); deltj++)
          {
       //       cout << feat << " "<< delti << " " << deltj << " obs " << _observed[feat][delti][deltj] << " delta  " <<_FM->getFeatureArraydelta(feat,delti,deltj,i,i) << endl;
            if(_gis)
            {
              if(_FM->getFeatureArraydelta(feat,delti,deltj,i,i)==1)
              {
                _observed[feat][delti][deltj]++;
              }
              else{

              }
            }
         //   else
        //    {
         //     if(_IM->getFeatureArraydelta(feat,delti,deltj,i,i)==1)
        //      {
          //      _observed[feat][delti][deltj]++;
       //       }
       //     }

          }

        }

      }

  }
  return _observed;
}

void IT::setFeatureArraylambda(int feati, int ilambdaX, int ilambdaY,double valuelambda){
  assert(feati<_systX.size());
  assert(_gis == true);
  _FM->setFeatureArraylambda(feati,ilambdaX,ilambdaY,valuelambda);
}


