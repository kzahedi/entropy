#ifndef __IS_PARAMETER_H__
#define __IS_PARAMETER_H__

namespace entropy
{
  namespace iterativescaling
  {
    class IsParameter
    {
      public:
        IsParameter()
        {
          lambdavalue    = 0.0;
          lambdadeltaval = 0.0;
          maxit          = -1;
          konv           = 0.0;
          konvtime       = false;
          test           = false;
          time           = false;
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
        bool       konvtime;
        int        seconds;
    };
  }
}

#endif // __IS_PARAMETER_H__
