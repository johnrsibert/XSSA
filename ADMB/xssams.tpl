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
  //const double LOG_TWO_M_PI = log(TWO_M_PI);
  const double LOG_M_PI = log(M_PI);
  
  // for finite differences, run
  // xssams -noinit -est -nr 10 -l2 10000000  -l3 10000000 &> xssams.out&
  // to avoid writing buffers; otherwise, run
  // xssams -noinit -est -nr 10 &> xssams.out&
  // seems to work without writing buffers

  int fexists(const adstring& filename)
  {
    std::ifstream ifile(filename);
    if (!ifile)
       return 0;
    else
       return 1;
  }


  template <typename SCALAR> SCALAR logit(const SCALAR& p)
  {
     SCALAR a = log(p/(1.0-p));
     return a;
  }
  template double logit<double>(const double& p);
  template dvariable logit<dvariable>(const dvariable& p);
  //template df1b2variable logit(const df1b2variable& p); 

  template <typename SCALAR> SCALAR alogit(const SCALAR& a)
  {
     SCALAR p = 1.0/(1.0+(mfexp(-a))+1e-20);
     return p;
  }
  template double alogit<double>(const double& a);
  template dvariable alogit<dvariable>(const dvariable& a);
  //template df1b2variable alogit<df1b2variable>(const df1b2variable& a); 
  /*
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

  init_int phase_T12;
  init_number init_T12;
  !! TTRACE(init_T12,phase_T12)

  init_int phase_T21;
  init_number init_T21;
  !! TTRACE(init_T21,phase_T21)

  init_int phase_Fmsy;
  init_number init_Fmsy;
  !! TTRACE(init_Fmsy,phase_Fmsy)

  init_int use_Fmsy_prior;
  init_number Fmsy_prior;
  init_number sdFmsy_prior;
  number logFmsy_prior;
  !! logFmsy_prior = log(Fmsy_prior);
  !! TRACE(use_Fmsy_prior)
  !! TTRACE(Fmsy_prior,sdFmsy_prior)
  number varFmsy_prior;
  !! varFmsy_prior = square(sdFmsy_prior);

  init_int phase_MSY;
  init_number init_MSY;
  !! TTRACE(init_MSY,phase_MSY)

  init_int phase_sdlogProc;
  init_number init_sdlogProc;
  !! TTRACE(init_sdlogProc,phase_sdlogProc)

  init_int phase_sdlogYield;
  init_number init_sdlogYield;
  !! TTRACE(init_sdlogYield,phase_sdlogYield)

  init_int phase_meanProportion_local;
  init_number init_meanProportion_local;
  !! TTRACE(init_meanProportion_local,phase_meanProportion_local)

  init_int phase_sdProportion_local;
  init_number init_sdProportion_local;
  !! TTRACE(init_sdProportion_local,phase_sdProportion_local)

  init_int phase_qProp;
  init_number init_qProp;
  !! TTRACE(init_qProp,phase_qProp)

  init_int use_robustY;
  init_int phase_pcon;
  init_number init_pcon;
  !! TTRACE(init_pcon,phase_pcon)

  int pininit;
  int lengthU;
  int utPop1;
  int utPop2;
  int userfun_entries;
  int status_blocks;
  ivector Fndxl;
  ivector Fndxu;
  int trace_init_pars;
  matrix residuals;
  ivector UU;
  ivector first_year;
  ivector last_year;

 LOCAL_CALCS;
    TRACE(ntime);
    TRACE(ngear);
    lengthU = ntime*(ngear+2);
    TRACE(lengthU)
    UU.allocate(1,lengthU);
    TTRACE(UU.indexmin(),UU.indexmax())
    for (int i = 1; i <= lengthU; i++)
       UU(i) = i;
    residuals.allocate(1,ntime,1,2*ngear+5);
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
    }
    clogf << "Zero catch bridging instances: " << nzero << endl;

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
    utPop1 = ngear*ntime;
    utPop2 = utPop1 + ntime;
    TTRACE(utPop1,utPop2)

    //if (1) ad_exit(1);

PARAMETER_SECTION
  // transfer parameters
  init_number logT12(phase_T12);
  init_number logT21(phase_T21);

  // logistic parameters
  //init_number logr(phase_r);
  //init_number logK(phase_K);

  init_number logFmsy(phase_Fmsy);
  init_number logMSY(phase_MSY);

  // general process error standard deviations
  init_number logsdlogProc(phase_sdlogProc);

  // observation error standard deviations
  init_number logsdlogYield(phase_sdlogYield);

  // logit transformed porportion local prior
  init_number LmeanProportion_local(phase_meanProportion_local);
  init_number logsdLProportion_local(phase_sdProportion_local);

  // non-linear term apportionment
  init_bounded_number qProp(0.0,1.0,phase_qProp);

  // robust yield likelihood proprtion contamination
  init_number Lpcon(phase_pcon);

  random_effects_vector U(1,lengthU);
  //vector U(1,lengthU);

  sdreport_number aT12;
  sdreport_number aT21;
  sdreport_number alogMSY;
  sdreport_number aMSY;
  sdreport_number aFmsy;
  sdreport_number ar;
  sdreport_number aK;
  sdreport_number asdlogProc;
  sdreport_number asdlogYield;
  sdreport_number aQ;

  objective_function_value nll;

PRELIMINARY_CALCS_SECTION
    userfun_entries = 0;
    status_blocks = 0;
    pininit = fexists(adstring(argv[0])+".pin");
    // set initial parameter value from data file
    TRACE(pininit)
    if (!pininit)
    {
       logT12 = log(init_T12+1.0e-20);
       logT21 = log(init_T21+1.0e-20);
       TTRACE(logT12,logT21)

       logFmsy = log(init_Fmsy);
       TTRACE(logFmsy,init_Fmsy)

       logMSY = log(init_MSY);
       TTRACE(logMSY,init_MSY)

       logsdlogYield = log(init_sdlogYield);
       logsdlogProc = log(init_sdlogProc);

       LmeanProportion_local = logit((double)init_meanProportion_local);
       logsdLProportion_local = log(logit((double)init_sdProportion_local));
       TTRACE(LmeanProportion_local,logsdLProportion_local)
       TTRACE(init_sdProportion_local,logit((double)init_sdProportion_local))
       //double prop = 1.0/(1.0+mfexp(-value(LmeanProportion_local)));
       double prop = alogit(value(LmeanProportion_local));
       TTRACE(LmeanProportion_local,prop)

       qProp = init_qProp;
       TTRACE(qProp,phase_qProp)

       if (!use_robustY)
       {
          phase_pcon = -1;
          init_pcon = 1e-25;
       }
       Lpcon = logit((const double&)init_pcon);

       double r = 2.0*mfexp(value(logFmsy));
       double K = 4.0*mfexp(value(logMSY))/(1.0e-20+r);
       //double K = immigrant_biomass[1];
       double Pop1 = prop*K;
       double Pop2 = K-Pop1;

       //dmatrix Ferr(1,ntime,1,ngear); Ferr.fill_randn(77);
       //dmatrix logPop1Err(1,ntime,1,2); logPop1Err.fill_randn(79);
       //dmatrix logPop2Err(1,ntime,1,2); logPop2Err.fill_randn(75);
       int ut = 0;
       TRACE(ut)
       for (int t = 1; t <= ntime; t++)
       {   
          for (int g = 1; g <= ngear; g++)
          {
             U(++ut) =  -5.0; //log(0.001);
          }
       }
       TTRACE(ut,utPop1)
       for (int t = 1; t <= ntime; t++)
       {   
          U(++ut) = log(Pop1);
       }
       TTRACE(ut,utPop2)
       for (int t = 1; t <= ntime; t++)
       {   
          U(++ut) = log(Pop2);
       }
       TRACE(ut)
       adstring pinname(adstring(argv[0])+".p00");
       ofstream pin(pinname);
       if (!pin)
       {
          cerr << "Error creating " << pinname << endl;
          ad_exit(1);
       }
       PINOUT(logT12)
       PINOUT(logT21)
       //PINOUT(logr)
       //PINOUT(logK)
       PINOUT(logFmsy)
       PINOUT(logMSY)
       PINOUT(logsdlogProc)
       PINOUT(logsdlogYield)
       PINOUT(LmeanProportion_local)
       PINOUT(logsdLProportion_local)
       PINOUT(alogit(value(Lpcon)))
       //PINOUT(U)
       pin << "# U:" << endl;
       pin << "#   F(t,g):" << endl;
       for (int tg = 1; tg <= (ntime*ngear); tg++)
           pin << " " << U(tg);
       pin << "\n#   Pop1(t):"<< endl;
       for (int tt = 1; tt <= ntime; tt++)
          pin << " " << U(utPop1+tt);
       pin << "\n#   Pop2(t):"<< endl;
       for (int tt = 1; tt <= ntime; tt++)
          pin << " " << U(utPop2+tt);
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
    TTRACE(utPop1,utPop2)
    TRACE(Fndxl)
    TRACE(Fndxu)
    TRACE(lengthU)
    TRACE(logT12)
    TRACE(logT21)
    TRACE(logFmsy)
    TRACE(logMSY)
    TRACE(logsdlogProc)
    TRACE(logsdlogYield)
    TRACE(logsdlogYield)
    TRACE(LmeanProportion_local)
    TRACE(logsdLProportion_local)
    TRACE(alogit(value(Lpcon)))
    TRACE(U)
    TTRACE(U(utPop1+1),U(utPop2+1))
    TRACE (trace_init_pars)
    clogf << endl;
    //if (1) ad_exit(1);

PROCEDURE_SECTION
  if (trace_init_pars)
  {
    clogf << "\nInitial step in to PROCEDURE_SECTION:" << endl;
    TRACE(userfun_entries)
    TTRACE(utPop1,utPop2)
    TRACE(Fndxl)
    TRACE(Fndxu)
    TRACE(lengthU)
    TRACE(logT12)
    TRACE(logT21)
    //TRACE(logr)
    //TRACE(logK)
    TRACE(logsdlogProc)
    TRACE(LmeanProportion_local)
    TRACE(logsdLProportion_local)
    TRACE(U)
    TTRACE(U(utPop1+1),U(utPop2+1))
    trace_init_pars = 0;
    TRACE (trace_init_pars)
    clogf << endl;
    //if (1) ad_exit(1);
  }

  nll = 0.0;

  step0(U(utPop1+1), U(utPop2+1), logsdlogProc, logFmsy, logMSY, LmeanProportion_local);

  for (int t = 2; t <= ntime; t++)
  {
     step(t, U(Fndxl(t-1),Fndxu(t-1)), U(Fndxl(t),Fndxu(t)), logsdlogProc,
             U(utPop1+t-1), U(utPop1+t), U(utPop2+t-1), U(utPop2+t),
             logFmsy, logMSY, logT12,logT21,LmeanProportion_local,logsdLProportion_local,qProp);
  }

  for (int t = 1; t <= ntime; t++)
  {
     obs(t,U(Fndxl(t),Fndxu(t)),U(utPop1+t-1),U(utPop1+t),
                                U(utPop2+t-1),U(utPop2+t),logsdlogYield,Lpcon);
  }

  if ((use_Fmsy_prior) && active(logFmsy) )
  {
     dvariable nll_Fmsy = 0.5*(log(TWO_M_PI*varFmsy_prior) + square(logFmsy - logFmsy_prior)/varFmsy_prior);
     nll += nll_Fmsy;
  }

  ar = 2.0*mfexp(logFmsy);
  aK = 4.0*mfexp(logMSY)/(1.0e-20+ar);
  aFmsy = mfexp(logFmsy);
  aMSY = mfexp(logMSY);
  //aBmsy = aK*(1.0-aFmsy/ar);
  asdlogProc = mfexp(logsdlogProc);
  aQ = qProp;
  aT12 = mfexp(logT12);
  aT21 = mfexp(logT21);
  //apcon = alogit((dvariable&)Lpcon);
  ++userfun_entries;
  int status_print = ntime;
  if (userfun_entries > lengthU)
     status_print = lengthU;
  if (userfun_entries % status_print == 0)
  //if (userfun_entries % lengthU == 0)
  {
     write_status(clogf);
  }


SEPARABLE_FUNCTION void step0(const dvariable& p11, const dvariable p21, const dvariable& lsdlogProc, const dvariable& lFmsy, const dvariable& lMSY, const dvariable& LmPropL) 
  // p11 U(utPop1+t-1) log N1 at start of time step
  // p21 U(utPop2+t-1) log N2 at start of time step

  // ensure that starting population size is near K
  dvariable r = 2.0*mfexp(lFmsy);
  dvariable K = 4.0*mfexp(lMSY)/(1.0e-20+r);
  dvariable PropL = alogit((dvariable&)LmPropL);
  //dvariable PropL = 1.0/(1.0+mfexp(-LmPropL));
  dvariable p10 = PropL*K;
  dvariable p20 = K-p10;
  //TTRACE(p10,p20)
  //TTRACE(log(p10),log(p20))
  //TTRACE(p11,p21)
  dvariable varlogPop = square(mfexp(lsdlogProc));
  dvariable Pnll = 0.0;
  Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(log(p10) - p11)/varlogPop);
  Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(log(p20) - p21)/varlogPop);

  nll += Pnll;


SEPARABLE_FUNCTION void step(const int t, const dvar_vector& f1, const dvar_vector& f2, const dvariable& lsdlogProc, const dvariable& p11, const dvariable p12, const dvariable& p21, const dvariable p22, const dvariable& lFmsy, const dvariable& lMSY, const dvariable& lT12, const dvariable& lT21, const dvariable& LmPropL, const dvariable& lsdLProportion_local, const dvariable& qP)
  // f1  U(Fndxl(t-1),Fndxu(t-1)) log F at start of time step
  // f2  U(Fndxl(t),Fndxu(t)      log F at end   of time step)
  // p11 U(utPop1+t-1) log N1 at start of time step
  // p12 U(utPop1+t)   log N1 at end   of time step
  // p21 U(utPop2+t-1) log N2 at start of time step
  // p22 U(utPop2+t)   log N2 at end   of time step


  dvariable varlogF = square(mfexp(lsdlogProc));
  dvariable varlogPop = square(mfexp(lsdlogProc));

  //dvariable r = mfexp(lr);
  //dvariable K = mfexp(lK);
  dvariable r = 2.0*mfexp(lFmsy);
  dvariable K = 4.0*mfexp(lMSY)/(1.0e-20+r);
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
     //TTRACE(g,Fnll)

     if (isnan(value(Fnll)) && (userfun_entries > 0))
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
  dvariable penN1 = 0.0;
  dvariable penN2 = 0.0;
  const double epsN = 1.0e-3;

  #define __FINITE_DIFFERENCE__
  #ifdef __FINITE_DIFFERENCE__
  #warning Building finite differnce approximation
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
  //TTRACE(nextLogN1,nextLogN2)
  dvariable prevN1;
  dvariable prevN2;
  dvariable dLogN1;
  dvariable dLogN2;
  double sdt = dt/niter;
  //TTRACE(niter,sdt)
  for (int ss = 1; ss <= niter; ss++)
  {
     prevN1 = mfexp(nextLogN1);
     prevN2 = mfexp(nextLogN2);
     dLogN1 = r*(1.0 - prevN1/K) - sumFg - T12 - 2.0*(1.0-q)*r*prevN2/K;
     dLogN2 = r*(1.0 - prevN2/K) - sumFg - T12 - 2.0*q*r*prevN1/K + T21*immigrant_biomass(t)/prevN2;

  ///dvariable nextLogN1 = log(posfun(nextN1,epsN,penN1));
     //nextLogN1 = nextLogN1 + dLogN1*sdt;
     //nextLogN2 = nextLogN2 + dLogN2*sdt;
     nextLogN1 = nextLogN1 + posfun(dLogN1,epsN,penN1)*sdt;
     nextLogN2 = nextLogN2 + posfun(dLogN2,epsN,penN2)*sdt;
     //TTRACE(nextLogN1,nextLogN2)
  }
  //TTRACE(mfexp(nextLogN1),mfexp(nextLogN2))
  //TTRACE(nextLogN1,nextLogN2)


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

  // #ifdef __FINITE_DIFFERENCE__
  #else 
  ///////////////////////////////////
  //
  // analytic integral
  // assumes dt = 1
  // Z evaluated at t-1
  //
  ///////////////////////////////////
  dvariable prevN1 = mfexp(p11+1e-8);
  dvariable prevN2 = mfexp(p21+1e-8);
  dvariable Z1 = sumFg + T12 + 2.0*(1.0-q)*(r/K)*prevN2;
  dvariable Z2 = sumFg + T12 + 2.0*q*(r/K)*prevN1;
  //TTRACE(Z1,Z2)
  //TTRACE(sumFg,T12)
  //TTRACE(prevN1,prevN2)
  // S[t,1] =        (K*(r-Z1))/(r+((K*(r-Z1)/S[t-1,1])-r)*exp(-(r-Z1)))
  dvariable nextN1 = (K*(r-Z1))/(r+((K*(r-Z1)/prevN1)-r)*mfexp(-(r-Z1)));
  // S[t,2] = (K*(r-Z2))/(r+((K*(r-Z2)/S[t-1,2])-r)*exp(-(r-Z2)*(1.0-T21)))
  dvariable num = (K*(r-Z2));
  dvariable den = (r+((K*(r-Z2)/prevN2)-r)*mfexp(-(r-Z2)*(1.0-T21*immigrant_biomass(t))));
  dvariable nextN2 = num/den;
  //dvariable nextN2 = (K*(r-Z2))/(r+((K*(r-Z2)/prevN2)-r)*mfexp(-(r-Z2)*(1.0-T21*immigrant_biomass(t))));
  dvariable nextLogN1 = log(posfun(nextN1,epsN,penN1));
  if (value(penN1) > 1.0e-8)
    TTRACE(nextN1,penN1)
  dvariable nextLogN2 = log(posfun(nextN2,epsN,penN2));
  if (value(penN2) > 1.0e-8)
    TTRACE(nextN2,penN2)
 
  if(userfun_entries>0)
  {
     if ( isnan(value(nextLogN1)) || isnan(value(nextLogN2)) ||
          isinf(value(nextLogN1)) || isinf(value(nextLogN2)) )
     {
        TTRACE(t,userfun_entries)
        TTRACE(prevN1,prevN2)
        TTRACE(r,K)
        TTRACE(sumFg,T12)
        TRACE(q)
        TTRACE(Z1,Z2)
        TTRACE(num,den)
        TRACE(K*(r-Z1))
        TRACE(r+((K*(r-Z2)/prevN2)-r)*mfexp(-(r-Z2)*(1.0-T21*immigrant_biomass(t))))
        TTRACE(nextN1,nextN2)
        TTRACE(nextLogN1,nextLogN2)
        write_status(clogf);
        ad_exit(1);
     }
  }
  #endif // #ifdef __FINITE_DIFFERENCE__

  dvariable Pnll = 0.0;
  Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(p12-nextLogN1)/varlogPop);
  Pnll += penN1;
  Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(p22-nextLogN2)/varlogPop);
  Pnll += penN2;

 
  // proportion local prior
  dvariable PLnll = 0.0;
  dvariable LpropL = nextLogN1 - nextLogN2;
  PLnll += 0.5*(log(TWO_M_PI*varLPropL) + square(LpropL - LmeanPropL)/varLPropL);

  if (isnan(value(PLnll)) && (userfun_entries>0))
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

  

SEPARABLE_FUNCTION void obs(const int t, const dvar_vector& f,const dvariable& pop11, const dvariable& pop21,  const dvariable& pop12, const dvariable& pop22, const dvariable& logsdlogYield, const dvariable& Lpc)
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

        dvariable pfat = alogit(Lpc);
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

        dvariable pcon = alogit(Lpc);
        dvariable b=2.*pcon/(width*sqrt(PI));  // This is the weight for the "robustifying" term

        dvariable norm_part = log(1.0-pcon)*mfexp(-diff2/(2.0*a2*v_hat));
   
        dvariable fat_part =  b/(1.+pow(diff2/(width2*a2*v_hat),2));
        Ynll += norm_part + fat_part;
     }

     else if (use_robustY == 3) // zero inflated normal
     {
        //dvariable pzero = alogit(Lpc(g));
        dvariable pzero = alogit((dvariable&)Lpc);
        if (obs_catch(t,g) > logZeroCatch)
        {
           Ynll += (1.0-pzero)*0.5*(log(TWO_M_PI*varlogYield) + square(obs_catch(t,g)-log_pred_yield(g))/varlogYield);
        }
        else
        {
           Ynll += pzero*0.5*(log(TWO_M_PI*varlogYield));
        }
     }

     else // default log-normal likelihood
     {
        Ynll += 0.5*(log(TWO_M_PI*varlogYield) + square(obs_catch(t,g)-log_pred_yield(g))/varlogYield);

        if (isnan(value(Ynll)) && (userfun_entries>0))
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
  double r = 2.0*mfexp(value(logFmsy));
  double K = 4.0*mfexp(value(logMSY))/(1.0e-20+r);
  residuals(t,++rc) = value(pop21);
  residuals(t,++rc) = value(pop22);
  residuals(t,++rc) = K;
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
    double r = 2.0*exp(value(logFmsy));
    double K = 4.0*exp(value(logMSY))/(1.0e-20+r);
    cout << "\n# Status block:" << status_blocks << endl;
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
    s << "# logFmsy = " << logFmsy << " (" << active(logFmsy) <<")" << endl;
    s << "#    Fmsy = " << mfexp(logFmsy) << endl;
    s << "#   Fmsy_prior = " << Fmsy_prior << " (" << (use_Fmsy_prior>0) << ")" << endl;
    s << "# sdFmsy_prior = " << sdFmsy_prior << endl;
    //s << "# logr = " << logr << " (" << active(logr) <<")" << endl;
    s << "#    r = " << r << endl;
    s << "# logMSY = " << logMSY << " (" << active(logMSY) <<")" << endl;
    s << "#    MSY = " << mfexp(logMSY) << endl;
    //s << "#    Bmsy = " << aBmsy << endl;
    //s << "# logK = " << logK << " (" << active(logK) <<")" << endl;
    s << "#    K = " << K << endl;
    //s << "#     logsdlogF: " << logsdlogF 
    //         <<  " (" << active(logsdlogF) <<")" << endl;
    //s << "#        sdlogF: " << mfexp(logsdlogF) << endl;
    s << "#   logsdlogProc: " << logsdlogProc
             <<  " (" << active(logsdlogProc) <<")" << endl;
    s << "#      sdlogPop: " << mfexp(logsdlogProc) << endl;
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
    s << "# pcon = " << alogit(value(Lpcon)) << " (" << active(Lpcon) <<")" << endl;
    s << "# qProp = " << qProp << " (" << active(qProp) << ")" << endl;
    // to keep the diagnostics R script happy
    s << "#     logsdlogF: " << logsdlogProc 
	     <<  " (" << active(logsdlogProc) <<")" << endl;
    s << "#        sdlogF: " << mfexp(logsdlogProc) << endl;
    s << "#   logsdlogPop: " << logsdlogProc
	     <<  " (" << active(logsdlogProc) <<")" << endl;
    s << "#      sdlogPop: " << mfexp(logsdlogProc) << endl;

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

