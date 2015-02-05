#include <TMB.hpp>
#include <iostream>
//#include "trace.h"
//std::ofstream clogf;//("xssams.log");
//clogf.open("xssams.log", ios::out);

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
   //std::cout << Fndxu << std::endl;
   //TRACE(Fndxu)

   PARAMETER(logT12);
   PARAMETER(logT21);
   PARAMETER(logr);
   PARAMETER(logK);
   PARAMETER_VECTOR(logsdlogF);
   PARAMETER(logsdlogPop);
   PARAMETER_VECTOR(logsdlogYield)
   PARAMETER(LmeanProportionLocal);
   PARAMETER(logsdLProportionLocal);
   PARAMETER_VECTOR(Lpfat);

   PARAMETER_VECTOR(U);
   Type r = exp(logr);
   Type K = exp(logK);
   vector <Type> varlogF = square(exp(logsdlogF));
   vector <Type>  sdlogF = exp(logsdlogF);
   Type varlogPop = square(exp(logsdlogPop));
   vector <Type> varlogYield = square(exp(logsdlogYield));
   Type T12 = exp(logT12);
   Type T21 = exp(logT21);

   array <Type> logF(ngear,ntime);
   for (int t = 0; t < maxtime; t++)
   {
       for (int g = 0; g < ngear; g++)
       {
          logF(g,t) = U(Fndxl(t)+g);
           
       }
   }


   Type nll = 0.0;
   for (int t = 1; t < maxtime; t++)
   {
      //step(t, U(Fndxl(t-1),Fndxu(t-1)), U(Fndxl(t),Fndxu(t)), logsdlogF,
      // p11 U(utPop1+t-1) log N1 at start of time step
      // p12 U(utPop1+t)   log N1 at end   of time step
      // p21 U(utPop2+t-1) log N2 at start of time step
      // p22 U(utPop2+t)   log N2 at end   of time step

      Type Fnll = 0.0;
      Type sumFg = 0.0;// total fishing mortality
      for (int g = 1; g <= ngear; g++)
      {
         Fnll += 0.5*(log(TWO_M_PI*varlogF) + square(logF(g,t)-logF(g,t-1))/varlogF);
         sumFg += exp(logF(g,t);
      }

      TYPE p21 = U(utPop2+t-1); // log N2 at start of time step
      TYPE p22 = U(utPop2+t);   // log N2 at start of time step
      Type nextLogN1 = U(utPop1+t-1); //p11;
      Type nextLogN2 = U(utPop1+t);   //p21;
      Type prevLogN1;
      Type prevLogN2;
      Type prevN1;
      Type prevN2;

      prevLogN1 = nextLogN1;
      prevLogN2 = nextLogN2;
      prevN1 = mfexp(prevLogN1);
      prevN2 = mfexp(prevLogN2);

      nextLogN1 += dt*(r*(1.0 - prevN1/K) - sumFg - T12 - r*prevN2/K);
      nextLogN2 += dt*(r*(1.0 - prevN2/K) - sumFg - T12 - r*prevN1/K + T21*immigrant_biomass(t)/prevN2);

      Type Pnll = 0.0;
      Type diff = log(mfexp(p12)+mfexp(p22))-log(mfexp(nextLogN1)+mfexp(nextLogN2));
      Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(diff)/varlogPop);


  // proportion local
  dvariable PLnll = 0.0;
  // Proportion_local = mfexp(nextLogN1)/(mfexp(nextLogN1)+mfexp(nextLogN2));
  // logit(p) = log(N1)-log(N2)
  dvariable LpropL = nextLogN1 - nextLogN2;
  //TTRACE(alogit(LmeanPropL),alogit(LpropL))
  PLnll += 0.5*(log(TWO_M_PI*varLPropL) + square(LpropL - LmeanPropL)/varLPropL);
  if (isnan(value(PLnll)))
  {
     TRACE(PLnll)
     TTRACE(nextLogN1,nextLogN2)
     TTRACE(prevLogN1,prevLogN2)
     TTRACE(prevN1,prevN2)
     TTRACE(r,K)
     TTRACE(LpropL,LmeanPropL)
     TRACE(varLPropL)
     write_status(clogf);
     ad_exit(1);
  }

  //clogf << t << " " << Fnll << " " <<Pnll << " " <<PLnll << " " << nll << endl;
  nll += (Fnll+Pnll+PLnll);


   }




   return nll;
}
