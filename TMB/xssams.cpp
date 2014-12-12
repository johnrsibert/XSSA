#include <TMB.hpp>

template <class Type> 
Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
   DATA_INTEGER(ntime);
   DATA_INTEGER(ngear);
   DATA_SCALAR(dt); // integration time step
   DATA_INTEGER(fr);
   //DATA_SCALAR(init_logT12);
   //DATA_INTEGER(phase_logT12);
   //DATA_SCALAR(init_logT21);
   //DATA_INTEGER(phase_logT21);
   //DATA_SCALAR(init_logr);
   //DATA_INTEGER(phase_logr);
   //DATA_SCALAR(init_logK);
   //DATA_INTEGER(phase_logK);
   //DATA_VECTOR(init_logsdlogF);
   //DATA_INTEGER(phase_logsdlogF);
   //DATA_VECTOR(init_logsdlogPop);
   //DATA_INTEGER(phase_logsdlogPop);
   //DATA_VECTOR(init_logsdlogYield);
   //DATA_INTEGER(phase_logsdlogYield);
   //DATA_SCALAR(init_LmeanProportion_local);
   //DATA_INTEGER(phase_LmeanProportion_local);
   //DATA_SCALAR(init_logsdLProportion_local);
   //DATA_INTEGER(phase_logsdLProportion_local);
   DATA_INTEGER(ss);      // 1/dt number of interations in the integration
   DATA_MATRIX(ObsCatch);
   DATA_VECTOR(ImmigrantBiomass)
   DATA_INTEGER(lengthU);
   DATA_INTEGER(utPop1);
   DATA_INTEGER(utPop2);
   DATA_INTEGER(maxtime);
   DATA_FACTOR(Fndxl);
   DATA_FACTOR(Fndxu);

   PARAMETER(logT12);
   PARAMETER(logT21);
   PARAMETER(logr);
   PARAMETER(logK);
   PARAMETER_VECTOR(logsdlogF);
   PARAMETER_VECTOR(logsdlogPop);
   PARAMETER_VECTOR(logsdlogYield)
   PARAMETER(LmeanProportionLocal);
   PARAMETER(logsdLProportionLocal);

   PARAMETER_VECTOR(U);
   Type r = exp(logr);
   Type K = exp(logK);
   vector <Type> varlogF = square(exp(logsdlogF));
   vector <Type>  sdlogF = exp(logsdlogF);
   vector <Type> varlogPop = square(exp(logsdlogPop));
   vector <Type> varlogYield = square(exp(logsdlogYield));
   Type T12 = exp(logT12);
   Type T21 = exp(logT21);

   Type nll = 0.0;
   for (int t = 0; t < ntime; t++)
   {
      for (int g = 0; g < ngear; g++)
      {
         nll += -dnorm(U(Fndxl(t-1)+g)-U(fndxl(t)+g),sdlogF);
      }

   }



   }

   return nll;
}
