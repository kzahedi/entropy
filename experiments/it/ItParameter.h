#ifndef __IT_PARAMETER_H__
#define __IT_PARAMETER_H__

class ItParameter
{
  public:
    ItParameter()
    {
      lambdavalue    = 0.0;
      lambdadeltaval = 0.0;
      maxit          = -1;
      konv           = 0.0;
      test           = false;
      time           = false;
      seconds        = 0;
      sigma          = 0;
    };
    // ~ItParameter();

    //ItParameter(const ItParameter);
    //ItParameter operator=(const ItParameter);

    double     lambdavalue;
    double     lambdadeltaval;
    double     sigma;
    int        maxit;
    double     konv;
    bool       test;
    bool       time;
    int        seconds;
};


#endif // __IT_PARAMETER_H__
