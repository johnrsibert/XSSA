#include <TMB.hpp>

//#undef REPORT
//#define REPORT(object) report << "# " << #object ":\n" << setprecision(5) << object << endl;

const double TWO_M_PI = 2.0 * M_PI;
const double LOG_TWO_M_PI = log(TWO_M_PI);
const double LOG_M_PI = log(M_PI);


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

template < class Type > Type xlogit(Type p)
{
   Type a = log(p / (1.0 - p));
   return a;
}

template < class Type > Type alogit(const Type & a)
{
   Type p = 1.0 / (1.0 + mfexp(-a));
   return p;
}


template < class Type > Type objective_function < Type >::operator()()
{
   DATA_INTEGER(ngear);
   DATA_INTEGER(ntime);
   DATA_SCALAR(dt);
   Type logdt = log(dt);
   DATA_MATRIX(obs_catch);
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
   REPORT(logsdlogYield);

   // logit transformed porportion local prior
   PARAMETER(LmeanProportion_local);
   REPORT(LmeanProportion_local);
   PARAMETER(logsdLProportion_local);
   REPORT(logsdLProportion_local);

   // non-linear term apportionment
   //init_bounded_number qProp(0.0,1.0,phase_qProp);
   PARAMETER(qProp);

   // robust yield likelihood proprtion contamination
   PARAMETER_VECTOR(Lpfat);	//(1,ngear,phase_pfat);

   //random_effects_vector U(1,lengthU);
   //vector U(1,lengthU);
   PARAMETER_VECTOR(U);

   matrix < Type > residuals(ntime, 2 * ngear + 5);

   //PRELIMINARY_CALCS_SECTION
   int userfun_entries = 0;
   int status_blocks = 0;

   //PROCEDURE_SECTION
   //objective_function_value nll;
   Type nll = 0.0;

   //  step0(U(utPop1+1), U(utPop2+1), logsdlogPop, logK, LmeanProportion_local);
   //SEPARABLE_FUNCTION void step0(const dvariable& p11, const dvariable p21, const dvariable& lsdlogPop, const dvariable& lK, const dvariable& LmPropL) 
   // p11 U(utPop1+t-1) log N1 at start of time step
   // p21 U(utPop2+t-1) log N2 at start of time step

   // ensure that starting population size is near K
   Type K = mfexp(logK);
   //dvariable PropL = alogit(LmPropL);
   //Type PropL = 1.0/(1.0+mfexp(-LmeanProportion_local));
   Type PropL = alogit(-LmeanProportion_local);
   Type p10 = PropL * K;
   Type p20 = K - p10;
   Type p11 = U(utPop1);
   Type p21 = U(utPop2);
   Type varlogPop = square(mfexp(logsdlogPop));
   Type Pnll = 0.0;
   Pnll +=
      0.5 * (log(TWO_M_PI * varlogPop) +
	     square(log(p10) - p11) / varlogPop);
   Pnll +=
      0.5 * (log(TWO_M_PI * varlogPop) +
	     square(log(p20) - p21) / varlogPop);

   nll += Pnll;


   //for (int t = 2; t <= ntime; t++)
   for (int t = 1; t < ntime; t++)
   {
      // step(t, U(Fndxl(t-1),Fndxu(t-1)), U(Fndxl(t),Fndxu(t)), logsdlogF,
      //         U(utPop1+t-1), U(utPop1+t), U(utPop2+t-1), U(utPop2+t), logsdlogPop,
      //         logr,logK,logT12,logT21,LmeanProportion_local,
      //         logsdLProportion_local,qProp);
      //SEPARABLE_FUNCTION void step(const int t, const dvar_vector& f1, const dvar_vector& f2, const dvariable& lsdlogF, const dvariable& p11, const dvariable p12, const dvariable& p21, const dvariable p22, const dvariable& lsdlogPop, const dvariable& lr, const dvariable& lK, const dvariable& lT12, const dvariable& lT21, const dvariable& LmPropL, const dvariable& lsdLProportion_local, const dvariable& qP)
      // f1  U(Fndxl(t-1),Fndxu(t-1)) log F at start of time step
      // f2  U(Fndxl(t),Fndxu(t)      log F at end   of time step)
      // p11 U(utPop1+t-1) log N1 at start of time step
      // p12 U(utPop1+t)   log N1 at end   of time step
      // p21 U(utPop2+t-1) log N2 at start of time step
      // p22 U(utPop2+t)   log N2 at end   of time step


      Type varlogF = square(mfexp(logsdlogF));
      Type varlogPop = square(mfexp(logsdlogPop));

      Type r = mfexp(logr);
      Type K = mfexp(logK);
      Type T12 = mfexp(logT12);
      Type T21 = mfexp(logT21);
      Type LmeanPropL = LmeanProportion_local;
      Type varLPropL = square(mfexp(logsdLProportion_local));

      Type sumFg = 0.0;		//sum(mfexp(ft1)); // total fishing mortality
      Type Fnll = 0.0;
      for (int g = 0; g < ngear; g++)
      {
	 // Fnll += 0.5*(log(TWO_M_PI*varlogF) + square(ft1(g)-ft2(g))/varlogF);
	 Fnll +=
	    0.5 * (log(TWO_M_PI * varlogF) +
		   square(U(Fndxl(t - 1) + g) -
			  U(Fndxl(t) + g)) / varlogF);
	 sumFg += mfexp(U(Fndxl(t - 1) + g));

	 if (isnan(Fnll))
	 {
	    //TTRACE(Fnll,varlogF)
	    //TTRACE(ft1(g),ft2(g))
	    //TTRACE(t,g)
	    //write_status(clogf);
	    //ad_exit(1);
	    return Fnll;
	 }
      }

      Type q = qProp;

      //       unstable for large r
      //Type prevN1 = mfexp(p11);
      //Type prevN2 = mfexp(p21);
      //Type nextLogN1 = p11 + dt*(r*(1.0 - prevN1/K) - sumFg - T12 - 2.0*(1.0-q)*r*prevN2/K);
      //Type nextLogN2 = p21 + dt*(r*(1.0 - prevN2/K) - sumFg - T12 - 2.0*q*r*prevN1/K + T21*immigrant_biomass(t)/prevN2);

      //       do multiple iterations per time step
      //       niter = 1 gives idential results to the code above
      int niter = 16;
      if (dt == 0.25)
      {
	 niter = 4;
      }
      //Type nextLogN1 = p11;
      Type nextLogN1 = U(utPop1 + t - 1);
      //Type nextLogN2 = p21;
      Type nextLogN2 = U(utPop2 + t - 1);
      Type prevN1;
      Type prevN2;
      Type dLogN1;
      Type dLogN2;
      Type sdt = dt / niter;
      for (int ss = 1; ss <= niter; ss++)
      {
	 prevN1 = mfexp(nextLogN1);
	 prevN2 = mfexp(nextLogN2);

	 dLogN1 =
	    r * (1.0 - prevN1 / K) - sumFg - T12 - 2.0 * (1.0 -
							  q) * r * prevN2 /
	    K;
	 dLogN2 =
	    r * (1.0 - prevN2 / K) - sumFg - T12 -
	    2.0 * q * r * prevN1 / K + T21 * immigrant_biomass(t) / prevN2;

	 nextLogN1 = nextLogN1 + dLogN1 * sdt;
	 nextLogN2 = nextLogN2 + dLogN2 * sdt;
      }

      /*
         if ( isnan(value(nextLogN1)) || isnan(value(nextLogN2)) ||
         isinf(value(nextLogN1)) || isinf(value(nextLogN2)) )
         {
         TTRACE(prevN1,prevN2)
         TTRACE(nextLogN1,nextLogN2)
         TTRACE(r,K)
         write_status(clogf);
         ad_exit(1);
         }
       */

      Type Pnll = 0.0;
      //Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(p12-nextLogN1)/varlogPop);
      Pnll +=
	 0.5 * (log(TWO_M_PI * varlogPop) +
		square(U(utPop1 + t) - nextLogN1) / varlogPop);
      //Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(p22-nextLogN2)/varlogPop);
      Pnll +=
	 0.5 * (log(TWO_M_PI * varlogPop) +
		square(U(utPop2 + t) - nextLogN2) / varlogPop);


      // proportion local prior
      Type PLnll = 0.0;
      Type LpropL = nextLogN1 - nextLogN2;
      PLnll +=
	 0.5 * (log(TWO_M_PI * varLPropL) +
		square(LpropL - LmeanPropL) / varLPropL);

      /*
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
       */

      //clogf << t << " " << Fnll << " " <<Pnll << " " <<PLnll << " " << nll << endl;
      nll += (Fnll + Pnll + PLnll);

   }				//for (int t = 1; t < ntime; t++)

   for (int t = 0; t < ntime; t++)
   {
      // obs(t,U(Fndxl(t),Fndxu(t)),U(utPop1+t-1),U(utPop1+t),
      //                            U(utPop2+t-1),U(utPop2+t),logsdlogYield,Lpfat);
      //SEPARABLE_FUNCTION void obs(const int t, const dvar_vector& f,const dvariable& pop11, const dvariable& pop21,  const dvariable& pop12, const dvariable& pop22, const dvariable& logsdlogYield, const dvar_vector& Lpf)
      // f  U(Fndxl(t),Fndxu(t)) 
      Type pop11 = U(utPop1 + t - 1);
      Type pop21 = U(utPop1 + t);
      Type pop12 = U(utPop2 + t - 1);
      Type pop22 = U(utPop2 + t);

      //dvar_vector ft(1,ngear); ///< log fishing mortality by gear
      //for (int g = 1; g <= ngear; g++)
      //{
      //   ft(g) = f(f.indexmin()+g-1);
      //}

      // sum of the average populations sizes over time step
      Type log_total_mean_pop;
      if (t < 2)
	 log_total_mean_pop = log(mfexp(pop21) +	// population 1
				  mfexp(pop22));	// population 2
      else
	 log_total_mean_pop = log(0.5 * (mfexp(pop11) + mfexp(pop21) +	// population 1
					 mfexp(pop12) + mfexp(pop22)));	// population 2

      //dvar_vector log_pred_yield(1,ngear);
      vector < Type > log_pred_yield(ngear);
      for (int g = 0; g < ngear; g++)
      {
	 log_pred_yield(g) = logdt + U(Fndxl(t) + g) + log_total_mean_pop;
      }

      Type varlogYield = square(mfexp(logsdlogYield));
      Type Ynll = 0.0;
      for (int g = 0; g < ngear; g++)
      {
	 // observation error
	 if (use_robustY == 1)	// normal + t distribution
	 {
	    Type z =
	       square(obs_catch(g, t) - log_pred_yield(g)) / varlogYield;

	    Type norm_part = 0.5 * (log(TWO_M_PI * varlogYield) + z);

	    // standard Cauchy density
	    Type fat_part = 1.0 / (M_PI * (1.0 + z));
	    Type pfat = alogit(Lpfat(g));
	    Ynll += log((1.0 - pfat) * mfexp(norm_part) + pfat * fat_part);
	 }
	 // this block is incorrect and should not be used
	 else if (use_robustY == 2)
	 {
	    //cerr << "Unsupported robustifation option: " << use_robustY << endl;
	    if (1)
	       return (1);
	    /*
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
	     */
	 } else			// default log-normal likelihood
	 {
	    Ynll +=
	       0.5 * (log(TWO_M_PI * varlogYield) +
		      square(obs_catch(g, t) -
			     log_pred_yield(g)) / varlogYield);

	    //if (isnan(value(Ynll)))
	    //{
	    //   TRACE(Ynll)
	    //   TTRACE(t,g)
	    //   write_status(clogf);
	    //   ad_exit(1);
	    //}
	 }
      }				// for (int g = 0; g < ngear; g++)
      nll += Ynll;
   REPORT(nll);

      /*
      ++userfun_entries;
      int status_print = ntime;
      if (userfun_entries > lengthU)
	 status_print = lengthU;
      //if (userfun_entries % status_print == 0)
      //{
      //   write_status(clogf);
      //}

      // dump stuff in "residual" matrix
      // matrix < Type > residuals(ntime, 2 * ngear + 5);
      int rc = -1;		// residuals column counter
      residuals(t, ++rc) = pop21;
      residuals(t, ++rc) = pop22;
      residuals(t, ++rc) = mfexp(logK);
      residuals(t, ++rc) = mfexp(logT21) * immigrant_biomass(t);
      Type propLa = pop21 - pop22;
      residuals(t, ++rc) = alogit(propLa);

      for (int g = 0; g < ngear; g++)
	 residuals(t, ++rc) = U(Fndxl(t) + g);
      for (int g = 0; g < ngear; g++)
	 residuals(t, ++rc) = log_pred_yield(g);
      //REPORT(residuals) 
      */
   } // for (int t = 0; t < ntime; t++)
   REPORT(nll);
   return nll;
}

/*
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
