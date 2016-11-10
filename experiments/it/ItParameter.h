#ifndef __IT_PARAMETER_H__
#define __IT_PARAMETER_H__

class IsParameter
{
  public:
    IsParameter()
    {
      lambdavalue    = 0.0;
      lambdadeltaval = 0.0;
      maxit          = -1;
      konv           = 0.0;
      test           = false;
      time           = false;
      konvtime       = false;
      seconds        = 0;
      sigma          = 0;
    };
    // ~IsParameter();

    //IsParameter(const IsParameter);
    //IsParameter operator=(const IsParameter);

    double     lambdavalue;
    double     lambdadeltaval;
    double     sigma;
    int        maxit;
    double     konv;
    bool       test;
    bool       time;
    bool 	   konvtime;
    int        seconds;
};


#endif // __IT_PARAMETER_H__
