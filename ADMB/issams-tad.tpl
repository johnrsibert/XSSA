GLOBALS_SECTION;
  #include <math.h>
  #include <adstring.hpp>
  #include "trace.h"
  #include <admodel.h>

  #undef PINOUT
  #define PINOUT(object) pin << "# " << #object ":\n" << setprecision(5) << object << endl;
  #undef REPORT
  #define REPORT(object) report << "# " << #object ":\n" << setprecision(5) << object << endl;

  ofstream clogf;
  const double TWO_M_PI = 2.0*M_PI;
  const double LOG_M_PI = log(M_PI);

   //  issams-tad -noinit -iprint 1 -shess &> issams.out&
   //  issams-tad -noinit -mcmc2 500000 -mcsave 20 -shess &> issams.out


  int fexists(const adstring& filename)
  {
    std::ifstream ifile(filename);
    if (!ifile)
       return 0;
    else
       return 1;
  }

  #include <df1b2fun.h>
  template <class Float>
  Float zilnorm(const double& pzero, Float var, const double& logobs, Float logpred)
  {
     const double TWO_M_PI = 2.0*M_PI;
     const double logZeroCatch = 0.0;
     Float ret;
     if (logobs > logZeroCatch)
        ret = (1.0-pzero)*0.5*(log(TWO_M_PI*var) + square(logobs-logpred)/var);
     else
        ret = pzero*0.5*(log(TWO_M_PI*var));
     return ret;
  }


TOP_OF_MAIN_SECTION
  arrmblsize = 50000000;
  gradient_structure::set_CMPDIF_BUFFER_SIZE(  150000000L);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(12550000L);
  gradient_structure::set_MAX_NVAR_OFFSET(3000000);

  adstring logname(adstring(argv[0])+"_program.log");
  clogf.open(logname);
  if ( !clogf ) {
    cerr << "Cannot open program log " << logname << endl;
    ad_exit(1);
  }
  cout << "Opened program log: " << logname << endl;
  pad();


DATA_SECTION
  !! streampos tmp = ad_comm::change_datafile_name("issams.dat");
  init_int nobs_gear;
  int ngear;
  init_int ntime;
  init_number dt;
  number logdt;
  !!  logdt = log(dt);
  init_matrix tcatch(1,nobs_gear,1,ntime);
  !!  TRACE(trans(tcatch))
  matrix obs_catch;

  init_int use_klingons;
  init_number klingon_multiplier;
  !! TTRACE(use_klingons,klingon_multiplier)
  !!  if (use_klingons)
  !!  {
  !!     ngear = nobs_gear + 1;
  !!     TTRACE(nobs_gear,ngear)
  !!     obs_catch.allocate(1,ntime,1,ngear);
  !!     for (int t = 1; t <= ntime; t++)
  !!     {
  !!        double tsum = 0.0;
  !!        for (int g = 1; g <= nobs_gear; g++) 
  !!        {
  !!           tsum += tcatch(g,t);
  !!           obs_catch(t,g) = tcatch(g,t);
  !!        }
  !!        obs_catch(t,ngear) = klingon_multiplier*tsum;
  !!    }
  !!  }
  !!  else
  !!  {
  !!     ngear = nobs_gear;
  !!     obs_catch = trans(tcatch);
  !!  }
  !!  TRACE(obs_catch)
  !! //if (1) exit(1);


  number logZeroCatch;

  init_matrix forcing_matrix(1,9,1,ntime);
  init_int fr;
  init_int use_mean_forcing;
  vector immigrant_biomass
  number mean_immigrant_biomass;
  number maximum_immigrant_biomass;
  !! immigrant_biomass = forcing_matrix(fr);
  !! mean_immigrant_biomass = mean(forcing_matrix(fr));
  !! maximum_immigrant_biomass = max(forcing_matrix(fr));
  !! TRACE(mean_immigrant_biomass)
  !! TRACE(maximum_immigrant_biomass)
  !! if (use_mean_forcing)
  !!    immigrant_biomass = mean_immigrant_biomass;
  !! TRACE(immigrant_biomass)

  init_int phase_Fmsy;
  init_number init_Fmsy;
  !! TTRACE(init_Fmsy,phase_Fmsy)

  init_int use_r_prior;
  init_number r_prior;
  init_number sdr_prior;
  number logr_prior;
  !! logr_prior = log(r_prior);
  !! TRACE(use_r_prior)
  !! TTRACE(r_prior,sdr_prior)
  number varr_prior;
  !! varr_prior = square(sdr_prior);
  !! TRACE(varr_prior)

  init_int phase_MSY;
  init_number init_MSY;
  !! TTRACE(init_MSY,phase_MSY)

  init_int phase_sdlogProc;
  init_number init_sdlogProc;
  !! TTRACE(init_sdlogProc,phase_sdlogProc)

  init_int phase_sdlogYield;
  init_number init_sdlogYield;
  !! TTRACE(init_sdlogYield,phase_sdlogYield)

  init_int use_Q;
  init_int phase_Q;
  init_number init_Q;
  !! TTRACE(init_Q,phase_Q)

  init_int use_robustY;
  !! TRACE(use_robustY)
  init_int phase_pcon;
  init_number pcon;

  int pininit;
  int lengthU;
  int utPop;
  int userfun_entries;
  int status_blocks;
  ivector Fndxl;
  ivector Fndxu;
  int trace_init_pars;
  matrix residuals;
  ivector UU;
  ivector first_year;
  ivector last_year;

  int nll_count;
  !! nll_count = 0;
  vector nll_vector(1,100000);
  !! nll_vector.initialize();

 LOCAL_CALCS;
    TRACE(ntime);
    TRACE(ngear);
    lengthU = ntime*(ngear+1);
    TRACE(lengthU)
    UU.allocate(1,lengthU);
    TTRACE(UU.indexmin(),UU.indexmax())
    for (int i = 1; i <= lengthU; i++)
       UU(i) = i;
    residuals.allocate(1,ntime,1,3*ngear+3);
    residuals.initialize();
    double ZeroCatch = 1.0; //1.0e-8;
    logZeroCatch = log(ZeroCatch);
    TTRACE(ZeroCatch,logZeroCatch)

    int nzero = ntime;
    if (use_robustY != 3)
    {
       int ziter = 0;
       while (nzero > 0)
       {
          ziter ++;
          nzero = 0;
          for (int g = 1; g <= ngear; g++)
            for (int t = 2; t <= ntime; t++)
               if ( (obs_catch(t,g) <= 0.0) 
                 && (obs_catch(t-1,g) > 0.0) && (obs_catch(t+1,g) > 0.0) )
               {
                  clogf << ++nzero << " " << ziter << endl;
                  obs_catch(t,g) = 0.5*(obs_catch(t-1,g) + obs_catch(t+1,g));
                  clogf << nzero<< " catch for gear " << g << " at time " << t
                     << " set to " << obs_catch(t,g)  << endl;
               }
       }
       clogf << "Zero catch bridging instances: " << nzero << endl;
    }

    first_year.allocate(1,ngear);
    first_year = ntime;
    last_year.allocate(1,ngear);
    last_year.initialize();
    ivector count_zero(1,ngear);
    count_zero.initialize();
    for (int t = 1; t<= ntime; t++)
    {
       for (int g = 1; g <= ngear; g++)
       {
          if ( (obs_catch(t,g) > 0.0) && (t < first_year(g)) )
             first_year(g) = t;
          if ( (obs_catch(t,g) > 0.0) && (t > last_year(g)) )
             last_year(g) = t;
       }
    }
    TRACE(first_year)
    TRACE(last_year)
    // dont't do this
    first_year = 1;
    last_year = ntime;
    TRACE(first_year)
    TRACE(last_year)

    for (int g = 1; g <= ngear; g++)
    {
       for (int t = first_year(g); t<= last_year(g); t++)
       {
          if (obs_catch(t,g) <= 0.0)
            count_zero(g) ++;
       }
    }
    TRACE(count_zero)

    dvector prop_zero(1,ngear);
    prop_zero = count_zero/double(ntime);
    double tprop_zero = (double)sum(count_zero)/(double)(ngear*ntime);
    TTRACE(prop_zero,tprop_zero)

    obs_catch = log(obs_catch+ZeroCatch);
    TRACE(obs_catch)
  
    // set up U indexing
    Fndxl.allocate(1,ntime);
    Fndxu.allocate(1,ntime);
    for (int t = 1; t<= ntime; t++)
    {
       Fndxl(t) = (t-1)*ngear+1;
       Fndxu(t) = t*ngear;
    }
    TRACE(Fndxl)
    TRACE(Fndxu)
    utPop = ngear*ntime;
    TRACE(utPop)

    //if (1) ad_exit(1);

PARAMETER_SECTION
  // logistic parameters
  init_number logFmsy(phase_Fmsy);
  init_number logMSY(phase_MSY);

  // general process error standard deviations
  init_number logsdlogProc(phase_sdlogProc);

  // observation error standard deviations
  init_number logsdlogYield(phase_sdlogYield);

  // abundance index proportionality constant
  init_number logQ(phase_Q);

  random_effects_vector U(1,lengthU);
  //vector U(1,lengthU);

  sdreport_number alogMSY;
  sdreport_number aMSY;
  sdreport_number alogFmsy;
  sdreport_number aFmsy;
  sdreport_number ar;
  sdreport_number aK;
  sdreport_number asdlogProc;
  sdreport_number asdlogYield;
  sdreport_number aQ;
  sdreport_number alogQ;

  objective_function_value nll;

PRELIMINARY_CALCS_SECTION
    userfun_entries = 0;
    status_blocks = 0;

    pininit = fexists(adstring(argv[0])+".pin");
    // set initial parameter value from data file
    TRACE(pininit)
    if (!pininit)
    {
       logFmsy = log(init_Fmsy);
       TTRACE(logFmsy,init_Fmsy)

       logMSY = log(init_MSY);
       TTRACE(logMSY,init_MSY)
       //logK = logMSY -logr + log(4.0);
       TTRACE(logFmsy,logMSY)

       logsdlogYield = log(init_sdlogYield);
       logsdlogProc = log(init_sdlogProc);
       logQ = log(init_Q);
       TTRACE(logQ,init_Q);

       double r = 2.0*mfexp(value(logFmsy));
       double K = 4.0*mfexp(value(logMSY))/(1.0e-20+r);

       int ut = 0;
       TRACE(ut)
       for (int t = 1; t <= ntime; t++)
       {   
          for (int g = 1; g <= ngear; g++)
          {
             U(++ut) = -5.0;
          }
       }
       TTRACE(ut,utPop)
       for (int t = 1; t <= ntime; t++)
       {   
          U(++ut) = log(K);
       }
       TRACE(ut)
       adstring pinname(adstring(argv[0])+".p00");
       ofstream pin(pinname);
       if (!pin)
       {
          cerr << "Error creating " << pinname << endl;
          ad_exit(1);
       }
       PINOUT(logFmsy)
       PINOUT(logMSY)
       PINOUT(logsdlogProc)
       PINOUT(logsdlogYield)
       PINOUT(logQ)
       //PINOUT(U)
       pin << "# U:" << endl;
       pin << "#   F(t,g):" << endl;
       for (int tg = 1; tg <= (ntime*ngear); tg++)
           pin << " " << U(tg);
       pin << "\n#   Pop(t):"<< endl;
       for (int tt = 1; tt <= ntime; tt++)
          pin << " " << U(utPop+tt);
       pin << endl;
       if (!pin)
       {
          cerr << "Error writing " << pinname << endl;
          ad_exit(1);
       }
       else
       {
          clogf << "Successfully created " << pinname << endl;
          cout << "Successfully created " << pinname << endl;
       }
    }
    trace_init_pars = 1;

    clogf << "\nAt end of PRELIMINARY_CALCS_SECTION:"<<endl;
    TRACE(userfun_entries)
    TRACE(utPop)
    TRACE(Fndxl)
    TRACE(Fndxu)
    TRACE(lengthU)


    TRACE(logFmsy)
    TRACE(logMSY)
    TRACE(logsdlogProc)
    TRACE(logsdlogYield)
    TRACE(logQ)
    TRACE(U)
    TRACE(U(utPop+1))
    TRACE (trace_init_pars)
    clogf << endl;
    //if (1) ad_exit(1);

PROCEDURE_SECTION

  if (trace_init_pars)
  {
    clogf << "\nInitial step in to PROCEDURE_SECTION:" << endl;
    TRACE(userfun_entries)
    TRACE(utPop)
    TRACE(Fndxl)
    TRACE(Fndxu)
    TRACE(lengthU)
    TRACE(logsdlogProc)
    TRACE(logsdlogYield)
    TRACE(U)
    TRACE(U(utPop+1))
    trace_init_pars = 0;
    TRACE (trace_init_pars)
    clogf << endl;
    //if (1) ad_exit(1);
  }

  nll = 0.0;
  nll_count = 0;

  step0(U(utPop+1), logsdlogProc, logFmsy, logMSY, logQ);
 
  for (int t = 2; t <= ntime; t++)
  {
     step(t, U(Fndxl(t-1),Fndxu(t-1)), U(Fndxl(t),Fndxu(t)),
             U(utPop+t-1), U(utPop+t), logsdlogProc,
             logFmsy, logMSY, logQ);
  }
  
  for (int t = 1; t <= ntime; t++)
  {
     obs(t,U(Fndxl(t),Fndxu(t)),U(utPop+t-1),U(utPop+t), logsdlogYield);
  }

  if (use_r_prior)
  {
     dvariable logr = log(2.0)+logFmsy;
     dvariable nll_r = 0.5*(log(TWO_M_PI*varr_prior) + square(logr - logr_prior)/varr_prior);
     nll += nll_r;
     NLL_TRACE(nll_r)
     NLL_TRACE(nll)
  }

  // compute sdreport_numbers
  ar = 2.0*mfexp(logFmsy);
  aK = 4.0*mfexp(logMSY)/(1.0e-20+ar);
  aFmsy = mfexp(logFmsy);
  alogFmsy = logFmsy;
  aMSY = mfexp(logMSY);
  alogMSY = logMSY;
  asdlogProc = mfexp(logsdlogProc);
  asdlogYield = mfexp(logsdlogYield);
  aQ = mfexp(logQ);
  alogQ = logQ;

  ++userfun_entries;
  if (!mc_phase())
  {
    int status_print = ntime;
    if (userfun_entries > lengthU)
       status_print = lengthU;
    if (userfun_entries % status_print == 0)
    {
       write_status(clogf);
    }
  }
  //write_nll_vector();


SEPARABLE_FUNCTION void step0(const dvariable& p11, const dvariable& lsdlogProc, const dvariable& lFmsy, const dvariable& lMSY, const dvariable& tlogQ)
  // p11 U(utPop+t-1) log N at start of time step

  // ensure that starting population size is near K
  dvariable r = 2.0*mfexp(lFmsy);
  dvariable K = 4.0*mfexp(lMSY)/(1.0e-20+r);
  dvariable p10 = log(K);
  dvariable lsdlogPop = lsdlogProc;
  dvariable varlogPop = square(mfexp(lsdlogPop));
  dvariable Pnll = 0.0;
  Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(p10 - p11)/varlogPop);
  NLL_TRACE(Pnll)
  NLL_TRACE(nll)
  nll += Pnll;

  dvariable Qnll = 0.0;
  if (use_Q)
  {
     dvariable varlogQ = square(mfexp(lsdlogProc));
     dvariable lnQ=tlogQ;
     dvariable logib = lnQ + log(immigrant_biomass(1));
     Qnll += 0.5*(log(TWO_M_PI*varlogQ) + square(logib-p11)/varlogQ);
     nll += Qnll;
  }
  NLL_TRACE(Qnll)
  NLL_TRACE(nll)



SEPARABLE_FUNCTION void step(const int t, const dvar_vector& f1, const dvar_vector& f2, const dvariable& p11, const dvariable p12, const dvariable& lsdlogProc, const dvariable& lFmsy, const dvariable& lMSY, const dvariable lnQ)
  // f1  U(Fndxl(t-1),Fndxu(t-1)) log F at start of time step
  // f2  U(Fndxl(t),Fndxu(t)      log F at end   of time step)
  // p11 U(utPop+t-1) log N at start of time step
  // p12 U(utPop+t)   log N at end   of time step


  dvariable lsdlogF = lsdlogProc;
  dvariable varlogF = square(mfexp(lsdlogF));
  dvariable lsdlogPop = lsdlogProc;
  dvariable varlogPop = square(mfexp(lsdlogPop));
 
  dvariable r = 2.0*exp(lFmsy);
  dvariable K = 4.0*exp(lMSY)/(1.0e-20+r);
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
  NLL_TRACE(Fnll)

  dvariable sumFg = sum(mfexp(ft1)); // total fishing mortality

  dvariable prevlogN = p11;
  dvariable prevN = mfexp(prevlogN);
  dvariable rmF = r - sumFg;
  dvariable ermF = mfexp(-1.0*rmF);
  dvariable Krmf = K*rmF;
  dvariable nextN = Krmf/(((Krmf/prevN)*ermF) - r*ermF +r); // 5
  dvariable nextLogN = log(nextN);

  if ( isnan(value(nextLogN)) )
  {
     TRACE(prevN)
     TRACE(nextLogN)
     TTRACE(r,K)
     write_status(clogf);
     ad_exit(1);
  }

  dvariable Pnll = 0.0;
  Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(p12-nextLogN)/varlogPop);
  NLL_TRACE(Pnll)

  nll += (Fnll+Pnll);
  NLL_TRACE(nll)

  dvariable Qnll = 0.0;
  if (use_Q)
  {
     dvariable varlogQ = square(mfexp(lsdlogProc));
     dvariable logib = lnQ + log(immigrant_biomass(t));
     Qnll += 0.5*(log(TWO_M_PI*varlogQ) + square(logib-nextLogN)/varlogQ);
     nll += Qnll;
  }
  NLL_TRACE(Qnll)
  NLL_TRACE(nll)

SEPARABLE_FUNCTION void obs(const int t, const dvar_vector& f,const dvariable& pop11, const dvariable& pop21, const dvariable& lsdlogYield)
  // f2  U(Fndxl(t),Fndxu(t)) 
  // pop11 U(utPop1+t-1)
  // pop21 U(utPop1+t)

  dvar_vector ft(1,ngear); ///< log fishing mortality by gear
  for (int g = 1; g <= ngear; g++)
  {
     ft(g) = f(f.indexmin()+g-1);
  }

  // sum of the average populations sizes over time step
  dvariable log_total_mean_pop;
  if (t < 2)
     log_total_mean_pop = pop21;
  else
     log_total_mean_pop = log( 0.5*(mfexp(pop11) + mfexp(pop21)) );

  dvar_vector log_pred_yield(1,ngear);
  for (int g = 1; g <= ngear; g++)
  {
     log_pred_yield(g) = logdt + ft(g) + log_total_mean_pop;
  }
  //TRACE(log_pred_yield)
  
  dvariable Ynll = 0.0;
  dvariable sdlogYield = mfexp(lsdlogYield);
  dvariable varlogYield = square(sdlogYield);
  for (int g = 1; g <= ngear; g++)
  {
     // observation error
     if (use_robustY == 1)  // normal + t distribution
     {
        dvariable z = square(obs_catch(t,g)-log_pred_yield(g))/varlogYield;

        dvariable norm_part = 0.5*(log(TWO_M_PI*varlogYield) + z);

        // standard Cauchy density
        dvariable fat_part = 1.0/(M_PI*(1.0 + z));
        Ynll += log((1.0-pcon)*mfexp(norm_part) + pcon*fat_part);
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

        dvariable b=2.*pcon/(width*sqrt(PI));  // This is the weight for the "robustifying" term

        dvariable norm_part = log(1.0-pcon)*mfexp(-diff2/(2.0*a2*v_hat));
   
        dvariable fat_part =  b/(1.+pow(diff2/(width2*a2*v_hat),2));
        Ynll += norm_part + fat_part;
     }

     else if (use_robustY == 3) // zero inflated normal
     {
        //TRACE(obs_catch(t,g))
        /*
        dvariable tmp;
        if (obs_catch(t,g) > logZeroCatch)
        {
           tmp = (1.0-pcon)*0.5*(log(TWO_M_PI*varlogYield) + square(obs_catch(t,g)-log_pred_yield(g))/varlogYield);
           //TTRACE(g,tmp)
           Ynll += tmp; //(1.0-pcon)*0.5*(log(TWO_M_PI*varlogYield) + square(obs_catch(t,g)-log_pred_yield(g))/varlogYield);
        }
        else
        {
           tmp = pcon*0.5*(log(TWO_M_PI*varlogYield));
           //TTRACE(g,tmp)
           Ynll += tmp; //pcon*0.5*(log(TWO_M_PI*varlogYield));
        }
        */
        Ynll += zilnorm((double&)pcon,varlogYield,obs_catch(t,g),(dvariable)log_pred_yield(g));
     }

     else // default normal likelihood
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
  } //for (int g = 1; g <= ngear; g++)

  nll += Ynll;
  NLL_TRACE(Ynll)
  NLL_TRACE(nll)

  // dump stuff in "residual" matrix
  int rc = 0; // residuals column counter
  double r = 2.0*mfexp(value(logFmsy));
  double K = 4.0*mfexp(value(logMSY))/(1.0e-20+r);
  residuals(t,++rc) = value(pop21);
  residuals(t,++rc) = K;
  residuals(t,++rc) = mfexp(value(logQ))*immigrant_biomass(t);

  for (int g = 1; g <= ngear; g++)
    residuals(t, ++rc) = value(ft(g));
  for (int g = 1; g <= ngear; g++)
    residuals(t, ++rc) = value(log_pred_yield(g));

FUNCTION void write_nll_vector(void)
    clogf << "\nnll_count " << nll_count << endl;
    for (int i = 1; i <= nll_count; i++)
       clogf << " " << setw(18) << setprecision(15) << nll_vector(i) 
                    << setw(6) << i << endl;

FUNCTION void write_status(ofstream& s)
    status_blocks ++;
    double r = 2.0*exp(value(logFmsy));
    double K = 4.0*exp(value(logMSY))/(1.0e-20+r);
    cout << "\n# Status block:" << status_blocks << endl;
    s << "\n# Status after "<< userfun_entries << " PROCEDURE_SECTION entries;" << endl;
    s << "# Status block " << status_blocks << endl;
    s << "# current phase = " << current_phase() << endl;
    s << "# nll = " << value(nll) << endl;
    s << "# nvar = " << initial_params::nvarcalc() << endl;
    s << "# logFmsy = " << logFmsy << " (" << active(logFmsy) <<")" << endl;
    s << "#    Fmsy = " << mfexp(logFmsy) << endl;
    s << "#    r = " << r << endl;
    s << "#   r_prior = " << r_prior << " (" << (use_r_prior>0) << ")" << endl;
    s << "# sdr_prior = " << sdr_prior << endl;
    s << "# logMSY = " << logMSY << " (" << active(logMSY) <<")" << endl;
    s << "#    MSY = " << mfexp(logMSY) << endl;
    s << "#    K = " << K << endl;
    s << "#   logsdlogProc: " << logsdlogProc
             <<  " (" << active(logsdlogProc) <<")" << endl;
    s << "#      sdlogProc: " << mfexp(logsdlogProc) << endl;
    s << "# logsdlogYield: " << logsdlogYield
             <<  " (" << active(logsdlogYield) <<")" << endl;
    s << "#    sdlogYield: " << mfexp(logsdlogYield) << endl;
    s << "#          logQ: " << logQ <<  " (" << active(logQ) <<")" << endl;
    s << "#             Q: " << mfexp(logQ) << " (" << use_Q <<")" << endl;
    s << "# pcon = " << pcon << endl;
    s << "# klingon_multiplier = " << klingon_multiplier << endl;
    // to keep the diagnostics R script happy
    s << "#     logsdlogF: " << logsdlogProc 
             <<  " (" << active(logsdlogProc) <<")" << endl;
    s << "#        sdlogF: " << mfexp(logsdlogProc) << endl;
    s << "#   logsdlogPop: " << logsdlogProc
             <<  " (" << active(logsdlogProc) <<")" << endl;
    s << "#      sdlogPop: " << mfexp(logsdlogProc) << endl;

    s << "# Residuals:" << endl;
    s << "  t    pop   K  forcing";
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
              << " " << residuals(t,3);
       int rc = 3;
       for (int g = 1; g <= 2*ngear; g++)
          s << " " << residuals(t,++rc);
       for (int g = 1; g <= ngear; g++)
          s << " " << obs_catch(t,g);
       s << endl;
    }



REPORT_SECTION
    clogf << "\nnll_count " << nll_count << endl;
    for (int i = 1; i <= nll_count; i++)
       clogf << " " << setw(18) << setprecision(15) << nll_vector(i) 
                    << setw(6) << i << endl;

    write_status(clogf);
    status_blocks --;
    report << "# Gmax = " << nll.gmax << endl;
    write_status(report);

