#include <TMB.hpp>
#include <iostream>
//#include <math.h>
//#include <adstring.hpp>
//#include "trace.h"
//#include <admodel.h>

//#undef PINOUT
//#define PINOUT(object) pin << "# " << #object ":\n" << setprecision(5) << object << endl;
//#undef REPORT
//#define REPORT(object) report << "# " << #object ":\n" << setprecision(5) << object << endl;

//ofstream clogf;
const double TWO_M_PI = 2.0*M_PI;
const double LOG_TWO_M_PI = log(TWO_M_PI);
const double LOG_M_PI = log(M_PI);


template <class Type> 
Type square(Type x)
{
   return x*x;
}

   /* 
  template <typename SCALAR> SCALAR logit(const SCALAR& p)
  {
     SCALAR a = log(p/(1.0-p));
     return a;
  }
  template double logit<double>(const double& p);
  //template dvariable logit<dvariable>(const dvariable& p);
  //template df1b2variable logit(const df1b2variable& p); 

  template <typename SCALAR> SCALAR alogit(const SCALAR& a)
  {
     SCALAR p = 1.0/(1.0+mfexp(-a));
     return p;
  }
  template double alogit<double>(const double& a);
  //template dvariable alogit<dvariable>(const dvariable& a);
  //template df1b2variable alogit(const df1b2variable& a); 

  dvector alogit(const dvector& a)
  {
     int n1 = a.indexmin();
     int n2 = a.indexmax();
     dvector p(n1,n2);
     for (int i = n1; i <= n2; i++)
        p(i) = alogit(a(i));
     return p;
  }
  */

template<class Type>
Type objective_function<Type>::operator() ()
{
   DATA_INTEGER(ngear);
   DATA_INTEGER(ntime);
   DATA_SCALAR(dt);
   DATA_MATRIX(ObsCatch);
   DATA_INTEGER(fr);
   DATA_VECTOR(immigrant_biomass);
   DATA_INTEGER(use_mean_forcing);
   DATA_INTEGER(phase_T12);
   DATA_SCALAR(init_T12);
   DATA_INTEGER(phase_T21);                 
   DATA_SCALAR(init_T21);
   DATA_INTEGER(phase_r);
   DATA_SCALAR(init_r);
   DATA_INTEGER(phase_K);
   DATA_SCALAR(init_K);
   DATA_INTEGER(phase_sdlogF);
   DATA_SCALAR(init_sdlogF);
   DATA_INTEGER(phase_sdlogPop);
   DATA_SCALAR(init_sdlogPop);
   DATA_INTEGER(phase_sdlogYield);
   DATA_SCALAR(init_sdlogYield);
   DATA_INTEGER(phase_meanProportion_local);
   DATA_SCALAR(init_meanProportion_local);
   DATA_INTEGER(phase_sdProportion_local);
   DATA_SCALAR(init_sdProportion_local);
   DATA_INTEGER(phase_qProp);
   DATA_SCALAR(init_qProp);
   DATA_INTEGER(use_robustY);
   DATA_INTEGER(phase_pfat);
   DATA_VECTOR(init_pfat);
   DATA_INTEGER(maxtime);
   DATA_INTEGER(lengthU);
   DATA_SCALAR(ss);
   DATA_IVECTOR(Fndxl);
   DATA_IVECTOR(Fndxu);
   DATA_INTEGER(utPop1);
   DATA_INTEGER(utPop2);

//PARAMETER_SECTION
   // transfer parameters
   PARAMETER(logT12);
   PARAMETER(logT21);
 
   // logistic parameters
   PARAMETER(logr);
   PARAMETER(logK);
 
   // random walk standard deviations
   PARAMETER(logsdlogF);
   PARAMETER(logsdlogPop);
 
   // observation error standard deviations
   PARAMETER(logsdlogYield);
 
   // logit transformed porportion local prior
   PARAMETER(LmeanProportion_local);
   PARAMETER(logsdLProportion_local);
 
   // non-linear term apportionment
   //init_bounded_number qProp(0.0,1.0,phase_qProp);
   PARAMETER(qProp);
 
   // robust yield likelihood proprtion contamination
   PARAMETER_VECTOR(Lpfat); //(1,ngear,phase_pfat);
 
   //random_effects_vector U(1,lengthU);
   //vector U(1,lengthU);
   PARAMETER_VECTOR(U);
 
   //objective_function_value nll;

//PRELIMINARY_CALCS_SECTION
   int userfun_entries = 0;
   int status_blocks = 0;

//PROCEDURE_SECTION
   Type nll = 0.0;

//  step0(U(utPop1+1), U(utPop2+1), logsdlogPop, logK, LmeanProportion_local);
//SEPARABLE_FUNCTION void step0(const dvariable& p11, const dvariable p21, const dvariable& lsdlogPop, const dvariable& lK, const dvariable& LmPropL) 
  // p11 U(utPop1+t-1) log N1 at start of time step
  // p21 U(utPop2+t-1) log N2 at start of time step

  // ensure that starting population size is near K
  Type K = exp(logK);
  //dvariable PropL = alogit(LmPropL);
  Type PropL = 1.0/(1.0+exp(-LmeanProportion_local));
  Type p10 = PropL*K;
  Type p20 = K-p10;
  Type p11 = U(utPop1);
  Type p21 = U(utPop2);
  Type varlogPop = square(exp(logsdlogPop));
  Type Pnll = 0.0;
  Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(log(p10) - p11)/varlogPop);
  Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(log(p20) - p21)/varlogPop);

  nll += Pnll;


  for (int t = 2; t <= ntime; t++)
  {
     step(t, U(Fndxl(t-1),Fndxu(t-1)), U(Fndxl(t),Fndxu(t)), logsdlogF,
             U(utPop1+t-1), U(utPop1+t), U(utPop2+t-1), U(utPop2+t), logsdlogPop,
             logr,logK,logT12,logT21,LmeanProportion_local,logsdLProportion_local,qProp);
  }

  for (int t = 1; t <= ntime; t++)
  {
     obs(t,U(Fndxl(t),Fndxu(t)),U(utPop1+t-1),U(utPop1+t),
                                U(utPop2+t-1),U(utPop2+t),logsdlogYield,Lpfat);
  }

  ++userfun_entries;
  int status_print = ntime;
  if (userfun_entries > lengthU)
     status_print = lengthU;
  //if (userfun_entries % status_print == 0)
  //{
  //   write_status(clogf);
  //}

/*


SEPARABLE_FUNCTION void step(const int t, const dvar_vector& f1, const dvar_vector& f2, const dvariable& lsdlogF, const dvariable& p11, const dvariable p12, const dvariable& p21, const dvariable p22, const dvariable& lsdlogPop, const dvariable& lr, const dvariable& lK, const dvariable& lT12, const dvariable& lT21, const dvariable& LmPropL, const dvariable& lsdLProportion_local, const dvariable& qP)
  // f1  U(Fndxl(t-1),Fndxu(t-1)) log F at start of time step
  // f2  U(Fndxl(t),Fndxu(t)      log F at end   of time step)
  // p11 U(utPop1+t-1) log N1 at start of time step
  // p12 U(utPop1+t)   log N1 at end   of time step
  // p21 U(utPop2+t-1) log N2 at start of time step
  // p22 U(utPop2+t)   log N2 at end   of time step


  dvariable varlogF = square(mfexp(lsdlogF));
  dvariable varlogPop = square(mfexp(lsdlogPop));

  dvariable r = mfexp(lr);
  dvariable K = mfexp(lK);
  dvariable T12 =mfexp(lT12);
  dvariable T21 =mfexp(lT21);
  dvariable LmeanPropL = LmPropL;
  dvariable varLPropL = square(mfexp(lsdLProportion_local));
  dvar_vector ft1(1,ngear);
  dvar_vector ft2(1,ngear);
  for (int g = 1; g <= ngear; g++)
  {
     ft1(g) = f1(f1.indexmin()+g-1);
     ft2(g) = f2(f2.indexmin()+g-1);
  }

  dvariable Fnll = 0.0;
  for (int g = 1; g <= ngear; g++)
  {
     Fnll += 0.5*(log(TWO_M_PI*varlogF) + square(ft1(g)-ft2(g))/varlogF);

     if (isnan(value(Fnll)))
     {
        TTRACE(Fnll,varlogF)
        TTRACE(ft1(g),ft2(g))
        TTRACE(t,g)
        write_status(clogf);
        ad_exit(1);
     }
  } 

  dvariable sumFg = sum(mfexp(ft1)); // total fishing mortality
  dvariable q = qP;

  //       unstable for large r
  //dvariable prevN1 = mfexp(p11);
  //dvariable prevN2 = mfexp(p21);
  //dvariable nextLogN1 = p11 + dt*(r*(1.0 - prevN1/K) - sumFg - T12 - 2.0*(1.0-q)*r*prevN2/K);
  //dvariable nextLogN2 = p21 + dt*(r*(1.0 - prevN2/K) - sumFg - T12 - 2.0*q*r*prevN1/K + T21*immigrant_biomass(t)/prevN2);

  //       do multiple iterations per time step
  //       niter = 1 gives idential results to the code above
  int niter = 16;
  if (dt == 0.25)
  {
     niter = 4;
  }
  dvariable nextLogN1 = p11;
  dvariable nextLogN2 = p21;
  dvariable prevN1;
  dvariable prevN2;
  dvariable dLogN1;
  dvariable dLogN2;
  double sdt = dt/niter;
  for (int ss = 1; ss <= niter; ss++)
  {
     prevN1 = mfexp(nextLogN1);
     prevN2 = mfexp(nextLogN2);

     dLogN1 = r*(1.0 - prevN1/K) - sumFg - T12 - 2.0*(1.0-q)*r*prevN2/K;
     dLogN2 = r*(1.0 - prevN2/K) - sumFg - T12 - 2.0*q*r*prevN1/K + T21*immigrant_biomass(t)/prevN2;

     nextLogN1 = nextLogN1 + dLogN1*sdt;
     nextLogN2 = nextLogN2 + dLogN2*sdt;
  }


  //       semi-implicit approximation
  //       this doesnt work because in the initial RE step
  //       most variables to zero includeing  p11 and p21
  //       so prevN1 and precN2 = 1
  //dvariable prevN1 = mfexp(p11);
  //dvariable prevN2 = mfexp(p21);
  //dvariable dtr = dt * r;
  //dvariable nextN1 =  prevN1*(1.0+dtr-dt*(sumFg+T12)-dtr*(1.0-q)*2.0*prevN2/K)/
  //                         (1.0+dtr*prevN1/K);
  //dvariable nextN2 = (prevN2*(1.0+dtr-dt*(sumFg+T12)-dtr*q*2.0*prevN1/K)+dt*T21*immigrant_biomass(t))/
  //                         (1.0+dtr*prevN2/K);

  if ( isnan(value(nextLogN1)) || isnan(value(nextLogN2)) ||
       isinf(value(nextLogN1)) || isinf(value(nextLogN2)) )
  {
     TTRACE(prevN1,prevN2)
     TTRACE(nextLogN1,nextLogN2)
     TTRACE(r,K)
     write_status(clogf);
     ad_exit(1);
  }

  dvariable Pnll = 0.0;
  Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(p12-nextLogN1)/varlogPop);
  Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(p22-nextLogN2)/varlogPop);

 
  // proportion local prior
  dvariable PLnll = 0.0;
  dvariable LpropL = nextLogN1 - nextLogN2;
  PLnll += 0.5*(log(TWO_M_PI*varLPropL) + square(LpropL - LmeanPropL)/varLPropL);

  if (isnan(value(PLnll)))
  {
     TRACE(PLnll)
     TTRACE(nextLogN1,nextLogN2)
     TTRACE(prevN1,prevN2)
     TTRACE(r,K)
     TTRACE(LpropL,LmeanPropL)
     TRACE(varLPropL)
     write_status(clogf);
     ad_exit(1);
  }

  //clogf << t << " " << Fnll << " " <<Pnll << " " <<PLnll << " " << nll << endl;
  nll += (Fnll+Pnll+PLnll);

  

SEPARABLE_FUNCTION void obs(const int t, const dvar_vector& f,const dvariable& pop11, const dvariable& pop21,  const dvariable& pop12, const dvariable& pop22, const dvariable& logsdlogYield, const dvar_vector& Lpf)
  // f2  U(Fndxl(t),Fndxu(t)) 
  // pop11 U(utPop1+t-1)
  // pop21 U(utPop1+t)
  // pop12 U(utPop2+t-1)
  // pop22 U(utPop2+t)

  dvar_vector ft(1,ngear); ///< log fishing mortality by gear
  for (int g = 1; g <= ngear; g++)
  {
     ft(g) = f(f.indexmin()+g-1);
  }

  // sum of the average populations sizes over time step
  dvariable log_total_mean_pop;
  if (t < 2)
     log_total_mean_pop = log(mfexp(pop21) +   // population 1
                              mfexp(pop22));   // population 2
  else
     log_total_mean_pop = log( 0.5*(mfexp(pop11) + mfexp(pop21) +   // population 1
                                    mfexp(pop12) + mfexp(pop22)) ); // population 2

  dvar_vector log_pred_yield(1,ngear);
  for (int g = 1; g <= ngear; g++)
  {
     log_pred_yield(g) = logdt + ft(g) + log_total_mean_pop;
  }

  dvariable varlogYield = square(mfexp(logsdlogYield));
  dvariable Ynll = 0.0;
  for (int g = 1; g <= ngear; g++)
  {
     // observation error
     if (use_robustY == 1)  // normal + t distribution
     {
        dvariable z = square(obs_catch(t,g)-log_pred_yield(g))/varlogYield;

        dvariable norm_part = 0.5*(log(TWO_M_PI*varlogYield) + z);

        // standard Cauchy density
        dvariable fat_part = 1.0/(M_PI*(1.0 + z));
        dvariable pfat = alogit(Lpf(g));
        Ynll += log((1.0-pfat)*mfexp(norm_part) + pfat*fat_part);
     }

     // this block is incorrect and should not be used
     else if (use_robustY == 2) 
     {
        cerr << "Unsupported robustifation option: " << use_robustY << endl;
        if (1)
           ad_exit(1);
        // based on a misreading of newreg2.cpp 
        double width=3.0;
        double width2=width*width;
        const double alpha = 0.7;
        const double a2 = square(alpha);
        dvariable diff2 = square(obs_catch(t,g)-log_pred_yield(g));
        dvariable v_hat = diff2+1.0e-80;

        dvariable pcon = alogit(Lpf(g));
        dvariable b=2.*pcon/(width*sqrt(PI));  // This is the weight for the "robustifying" term

        dvariable norm_part = log(1.0-pcon)*mfexp(-diff2/(2.0*a2*v_hat));
   
        dvariable fat_part =  b/(1.+pow(diff2/(width2*a2*v_hat),2));
        Ynll += norm_part + fat_part;
     }
     else // default log-normal likelihood
     {
        Ynll += 0.5*(log(TWO_M_PI*varlogYield) + square(obs_catch(t,g)-log_pred_yield(g))/varlogYield);

        if (isnan(value(Ynll)))
        {
           TRACE(Ynll)
           TTRACE(t,g)
           write_status(clogf);
           ad_exit(1);
        }
     }
  }

  nll += Ynll;

  // dump stuff in "residual" matrix
  int rc = 0; // residuals column counter
  residuals(t,++rc) = value(pop21);
  residuals(t,++rc) = value(pop22);
  residuals(t,++rc) = mfexp(value(logK));
  residuals(t,++rc) = mfexp(value(logT21))*immigrant_biomass(t);
  double propLa = value(pop21) - value(pop22);
  residuals(t,++rc) = alogit(propLa);

  for (int g = 1; g <= ngear; g++)
    residuals(t, ++rc) = value(ft(g));
  for (int g = 1; g <= ngear; g++)
    residuals(t, ++rc) = value(log_pred_yield(g));

FUNCTION void write_status(ofstream& s)
    double prop = alogit(value(LmeanProportion_local));
    status_blocks ++;
    cout << "\n# Status block " << status_blocks << endl;
    s << "\n# Status after "<< userfun_entries << " PROCEDURE_SECTION entries;" << endl;
    s << "# Status block " << status_blocks << endl;
    s << "# nll = " << value(nll) << endl;
    //s << "# maxG = " << nll.gmax << endl;
    s << "# nvar = " << initial_params::nvarcalc() << endl;
    s << "# current phase = " << current_phase() << endl;
    s << "#  logT12 = " << logT12 << " (" << active(logT12) <<")" << endl;
    s << "#  T12 = " << mfexp(logT12) << endl;
    s << "#  logT21 = " << logT21 << " (" << active(logT21) <<")" << endl;
    s << "#  T21 = " << mfexp(logT21) << endl;
    s << "# logr = " << logr << " (" << active(logr) <<")" << endl;
    s << "#    r = " << mfexp(logr) << endl;
    s << "# logK = " << logK << " (" << active(logK) <<")" << endl;
    s << "#    K = " << mfexp(logK) << endl;
    s << "#     logsdlogF: " << logsdlogF 
             <<  " (" << active(logsdlogF) <<")" << endl;
    s << "#        sdlogF: " << mfexp(logsdlogF) << endl;
    s << "#   logsdlogPop: " << logsdlogPop
             <<  " (" << active(logsdlogPop) <<")" << endl;
    s << "#      sdlogPop: " << mfexp(logsdlogPop) << endl;
    s << "# logsdlogYield: " << logsdlogYield
             <<  " (" << active(logsdlogYield) <<")" << endl;
    s << "#    sdlogYield: " << mfexp(logsdlogYield) << endl;
    s << "# LmeanProportion_local = " << LmeanProportion_local << " (" 
                               << active(LmeanProportion_local) <<")" << endl;
    s << "#                  prop = " << prop << endl;
    s << "# logsdLProportion_local = " << logsdLProportion_local<< " (" 
                                << active(logsdLProportion_local) <<")" << endl;
    s << "#    sdLProportion_local = " << mfexp(logsdLProportion_local) << endl;
    s << "#     sdProportion_local = " << alogit(value(mfexp(logsdLProportion_local))) << endl;
    s << "# pfat = " << alogit(value(Lpfat)) << " (" << active(Lpfat) <<")" << endl;
    s << "# qProp = " << qProp << " (" << active(qProp) << ")" << endl;
    s << "# Residuals:" << endl;
    s << "  t    pop1   pop2      K  forcing  propL";
    for (int g = 1; g <= ngear; g++)
       s << "     F" << g;
    for (int g = 1; g <= ngear; g++)
       s << "  predC" << g;
    for (int g = 1; g <= ngear; g++)
       s << "   obsC" << g;
    s << endl;
    for (int t = 1; t <= ntime; t++)
    {
       s << t << " " << residuals(t,1) << " " << residuals(t,2)
              << " " << residuals(t,3) << " " << residuals(t,4)
              << " " << residuals(t,5);
       int rc = 5;
       for (int g = 1; g <= 2*ngear; g++)
          s << " " << residuals(t,++rc);
       for (int g = 1; g <= ngear; g++)
          s << " " << obs_catch(t,g);
       s << endl;
    }



REPORT_SECTION
    write_status(clogf);
    status_blocks --;
    write_status(report);
*/
