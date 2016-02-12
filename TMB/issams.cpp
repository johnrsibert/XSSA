#include <TMB.hpp>
#include "trace.h"
#include <math.h>
const double TWO_M_PI = 2.0*M_PI;
const double LOG_M_PI = log(M_PI);
const double logZeroCatch = 1.0;

template < class Type > Type square(Type x)
{
   return x * x;
}

template < class Type > Type mfexp(Type x)
{
   double b = 60;
   if (x <= b && x >= -b)
   {
      return exp(x);
   } else if (x > b)
   {
      return exp(b) * (1. + 2. * (x - b)) / (1. + x - b);
   } else
   {
      return exp(-b) * (1. - x - b) / (1. + 2. * (-x - b));
   }
}
/*
template < class Type > Type logit(Type p)
{
   Type a = log(p / (Type(1.0) - p));
   return a;
}
*/

template < class Type > Type alogit(const Type & a)
{
   Type p = Type(1.0) / (Type(1.0) + mfexp(-a));
   return p;
}


//DATA_SECTION
template < class Type > Type objective_function < Type >::operator()()
{
   DATA_INTEGER(ngear);
   REPORT(ngear)
   DATA_INTEGER(ntime);
   REPORT(ntime);
   DATA_SCALAR(dt);
   REPORT(dt)
   Type logdt = log(dt);
   REPORT(logdt)
   DATA_MATRIX(obs_catch);
   REPORT(obs_catch);
   DATA_INTEGER(fr);
   REPORT(fr)
   DATA_VECTOR(immigrant_biomass);
   REPORT(immigrant_biomass);
   DATA_INTEGER(use_mean_forcing);
   REPORT(use_mean_forcing);
   DATA_INTEGER(phase_Fmsy);
   REPORT(phase_Fmsy);
   DATA_SCALAR(init_Fmsy);
   REPORT(init_Fmsy);
   DATA_INTEGER(use_r_prior);
   REPORT(use_r_prior);
   DATA_SCALAR(logr_prior);
   REPORT(logr_prior);
   DATA_SCALAR(varr_prior);
   REPORT(varr_prior);
   DATA_INTEGER(phase_MSY);
   REPORT(phase_MSY);
   DATA_SCALAR(init_MSY);
   REPORT(init_MSY);
   DATA_INTEGER(phase_sdlogProc);
   REPORT(phase_sdlogProc);
   DATA_SCALAR(init_sdlogProc);
   REPORT(init_sdlogProc);
   DATA_INTEGER(phase_sdlogYield);
   REPORT(phase_sdlogYield);
   DATA_SCALAR(init_sdlogYield);
   REPORT(init_sdlogYield);
   DATA_INTEGER(use_Q);
   REPORT(use_Q);
   DATA_SCALAR(init_Q);
   REPORT(init_Q);
   DATA_INTEGER(use_robustY);
   REPORT(use_robustY);
   DATA_SCALAR(pcon)
   REPORT(pcon)
   DATA_INTEGER(lengthU);
   REPORT(lengthU);
   DATA_IVECTOR(Fndxl);
   REPORT(Fndxl);
   DATA_IVECTOR(Fndxu);
   REPORT(Fndxu);
   DATA_INTEGER(utPop);
   REPORT(utPop);
   DATA_IVECTOR(first_year);
   REPORT(first_year);
   DATA_IVECTOR(last_year);
   REPORT(last_year)

   // logistic parameters
   PARAMETER(logFmsy);
   REPORT(logFmsy);
   PARAMETER(logMSY);
   REPORT(logMSY);

   // general process error standard deviations
   PARAMETER(logsdlogProc);
   REPORT(logsdlogProc);

   // observation error standard deviations
   PARAMETER(logsdlogYield);
   REPORT(logsdlogYield);

   // abundance index proportionality constant
   PARAMETER(logQ);
   REPORT(logQ);

   //random_effects_vector U(1,lengthU);
   PARAMETER_VECTOR(U);

   //residuals.allocate(1,ntime,1,3*ngear+3);
   matrix < double > residuals(ntime, 3 * ngear + 3);

   int userfun_entries = 0;
   int status_blocks = 0;

   Type nll = 0.0;

   // special case for initial population
   // ensure that starting population size is near K
   Type r = 2.0*mfexp(logFmsy);
   Type K = 4.0*mfexp(logMSY)/(1.0e-20+r);
   Type p10 = log(K);
   Type p11 = U(utPop);
   Type lsdlogPop = logsdlogProc;
   Type varlogPop = square(mfexp(lsdlogPop));
   Type Pnll = 0.0;
   Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(p10 - p11)/varlogPop);
   if (isnan(value(Pnll)))
       TRACE(Pnll)
   nll += Pnll;

   if (use_Q)
   {
      Type Qnll = 0.0;
      Type varlogQ = square(mfexp(logsdlogProc));
      Type logib = logQ + log(immigrant_biomass(0));
      Qnll += 0.5*(log(TWO_M_PI*varlogQ) + square(logib-p11)/varlogQ);
      if (isnan(value(Qnll)))
         TRACE(Qnll)
      nll += Qnll;
   }


   // update population until the end of time
   for (int t = 1; t < ntime; t++)
   {
      Type lsdlogF = logsdlogProc;
      Type varlogF = square(mfexp(lsdlogF));
      Type lsdlogPop = logsdlogProc;
      Type varlogPop = square(mfexp(lsdlogPop));
    
      Type r = 2.0*exp(logFmsy);
      Type K = 4.0*exp(logMSY)/(1.0e-20+r);
   
      Type Fnll = 0.0;
      Type sumFg = 0.0; // total fishing mortality
      for (int g = 0; g < ngear; g++)
      {
         Type ft1 = U(Fndxl(t-1)+g);
         Type ft2 = U(Fndxl(t)+g);
         sumFg += mfexp(ft1);
         Fnll += 0.5*(log(TWO_M_PI*varlogF) + square(ft1-ft2)/varlogF);
         if (isnan(value(Fnll)))
            TRACE(Fnll)
      } 
   
   
      Type p11 = U(utPop+t-1);
      Type p12 = U(utPop+t);
      Type prevLogN = p11;
      Type prevN = mfexp(prevLogN);
      Type rmF = r - sumFg;
      Type ermF = mfexp(-1.0*rmF);
      Type Krmf = K*rmF;
      Type nextN = Krmf/(((Krmf/prevN)*ermF) - r*ermF +r); // 5
      Type nextLogN = log(nextN);
   
      if ( isnan(value(nextLogN)) )
      {
          //write_status(clogf);
          TTRACE(prevLogN,nextLogN)
      }
   
      Type Pnll = 0.0;
      Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(p12-nextLogN)/varlogPop);
      if (isnan(value(Pnll)))
         TRACE(Pnll)
      nll += (Fnll+Pnll);
      if (use_Q)
      {
         Type Qnll = 0.0;
         Type varlogQ = square(mfexp(logsdlogProc));
         Type logib = logQ + log(immigrant_biomass(t));
         Qnll += 0.5*(log(TWO_M_PI*varlogQ) + square(logib-nextLogN)/varlogQ);
         if (isnan(value(Qnll)))
            TRACE(Qnll)
         nll += Qnll;
      }
   } //for (int t = 1; t < ntime; t++)

   // obs
   // loop through observations to compute likelihood
   for (int t = 0; t < ntime; t++)
   {
     Type pop11 = U(utPop+t-1);
     Type pop21 = U(utPop+t);

     vector <Type> ft(ngear);
     for (int g = 0; g < ngear; g++)
     {
        ft(g) = U(Fndxl(t)+g);
     }
     //TTRACE(t,ft)

     // sum of the average populations sizes over time step
     Type log_total_mean_pop;
     if (t < 2)
        log_total_mean_pop = pop21;
     else
        log_total_mean_pop = log( 0.5*(mfexp(pop11) + mfexp(pop21)) );

     vector <Type> log_pred_yield(ngear);
     for (int g = 0; g < ngear; g++)
     {
        log_pred_yield(g) = logdt + ft(g) + log_total_mean_pop;
     }
     //TTRACE(t,log_pred_yield)
   
     Type Ynll = 0.0;
     for (int g = 0; g < ngear; g++)
     {
        // observation error
        if ( (t >= first_year(g)) && (t <= last_year(g)) )
        {
           //log_pred_yield(g) = logdt + ft(g) + log_total_mean_pop;

           Type varlogYield = square(mfexp(logsdlogYield));
           if (use_robustY == 1)  // normal + t distribution
           {
              Type z = square(obs_catch(t,g)-log_pred_yield(g))/varlogYield;

              Type norm_part = 0.5*(log(TWO_M_PI*varlogYield) + z);

              // standard Cauchy density
              Type fat_part = 1.0/(M_PI*(1.0 + z));
              //Type pfat = alogit(Lpcon);
              Ynll += log((1.0-pcon)*mfexp(norm_part) + pcon*fat_part);
           }

           // this block is incorrect and should not be used
           else if (use_robustY == 2) 
           {
              std::cerr << "Unsupported robustifation option: " << use_robustY << std::endl;
              if (1)
                 exit(1);
           }

           else if (use_robustY == 3) // zero inflated normal
           {
              if (obs_catch(t,g) > logZeroCatch)
              {
                 Ynll += (1.0-pcon)*0.5*(log(TWO_M_PI*varlogYield) + square(obs_catch(t,g)-log_pred_yield(g))/varlogYield);
              } 
              else
              {
                 Ynll += pcon*0.5*(log(TWO_M_PI*varlogYield));
              }
              if (isnan(value(Ynll)))
                 TRACE(Ynll)
           }

           else // default normal likelihood
           {
              Ynll += 0.5*(log(TWO_M_PI*varlogYield) + square(obs_catch(t,g)-log_pred_yield(g))/varlogYield);

           }
         }
         else
         {
            log_pred_yield(g) = logZeroCatch;
            //HERE
         }
     }

     nll += Ynll;

     // dump stuff in "residual" matrix
     int rc = -1; // residuals column counter
     double r = 2.0*exp(value(logFmsy));
     double K = 4.0*exp(value(logMSY))/(1.0e-20+r);
     residuals(t,++rc) = value(pop21);
     residuals(t,++rc) = K;
     double Q = exp(value(logQ));
   //residuals(t,++rc) = Q; 
     residuals(t,++rc) = Q*immigrant_biomass(t); // produces error


     for (int g = 0; g < ngear; g++)
        residuals(t, ++rc) = value(ft(g));
     for (int g = 0; g < ngear; g++)
        residuals(t, ++rc) = value(log_pred_yield(g));

   } //for (int t = 0; t < ntime; t++)

   REPORT(residuals)

   if (use_r_prior)
   {
      Type logr = log(2.0)+logFmsy;
      Type nll_r = 0.5*(log(TWO_M_PI*varr_prior) + square(logr - logr_prior)/varr_prior);
      nll += nll_r;
      if (isnan(value(nll_r)))
         TRACE(nll_r)
   }

   return nll;
}

