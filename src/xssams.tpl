GLOBALS_SECTION;
  #include <math.h>
  #include <adstring.hpp>
  #include "trace.h"
  #include <fvar.hpp>
  #include <admodel.h>

  #undef PINOUT
  #define PINOUT(object) pin << "# " << #object ":\n" << setprecision(5) << object << endl;
  #undef REPORT
  #define REPORT(object) report << "# " << #object ":\n" << setprecision(5) << object << endl;

  ofstream clogf;
  const double TWO_M_PI = 2.0*M_PI;
  const double LOG_TWO_M_PI = log(TWO_M_PI);
  const double LOG_M_PI = log(M_PI);

  // ./xssams -noinit -nr 2 -est -l2 10000000  -l3 10000000
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


DATA_SECTION;
  init_adstring catch_dat_file;
  init_adstring forcing_dat_file;

  init_int ntime;
  init_int ngear;
  init_number dt; // integration time step
  init_int fr;
  !! TTRACE(dt,fr)

  init_number init_logT12;
  init_int phase_logT12;
  !! TTRACE(init_logT12,phase_logT12)


  init_number init_logT21;
  init_int phase_logT21;
  !! TTRACE(init_logT21,phase_logT21)

  init_number init_logr;
  init_int phase_logr;
  !! TTRACE(init_logr,phase_logr)

  init_number init_logK;
  init_int phase_logK;
  !! TTRACE(init_logK,phase_logK)

  init_number init_logsdlogF;
  init_int phase_logsdlogF;
  !! TTRACE(init_logsdlogF,phase_logsdlogF)

  init_vector init_logsdlogPop(1,2);
  init_int phase_logsdlogPop;
  !! TTRACE(init_logsdlogPop,phase_logsdlogPop)

  init_vector init_logsdlogYield(1,ngear);
  init_int phase_logsdlogYield;
  !! TTRACE(init_logsdlogYield,phase_logsdlogYield)

  init_number init_LmeanProportion_local;
  init_int phase_LmeanProportion_local;
  !! TTRACE(init_LmeanProportion_local,phase_LmeanProportion_local)

  init_number init_logsdLProportion_local;
  init_int phase_logsdLProportion_local;
  !! TTRACE(init_logsdLProportion_local,phase_logsdLProportion_local)

  init_int use_robustF;
  init_vector init_pfat(1,ngear);
  init_int phase_pfat;
  !! TTRACE(init_pfat,phase_pfat)

  number ss;      // 1/dt number of interations in the integration
  matrix obs_catch;
  matrix forcing_matrix;
  vector immigrant_biomass;
  number mean_immigrant_biomass;
  int pininit;
  number Fmsy;
  number Bmsy;
  number MSY;
  number ZeroCatch;
  number logZeroCatch;
  number ZeroF;
  number logZeroF;
  int lengthU;
  int utPop1;
  int utPop2;
  int userfun_entries;
  int maxtime;
  ivector Fndxl;
  ivector Fndxu;
  int trace_init_pars;
  matrix residuals;
  ivector UU;

 LOCAL_CALCS;
    maxtime = ntime;
    TTRACE(ntime,maxtime);
    TRACE(ngear);
    lengthU = ntime*(ngear+2);
    TRACE(lengthU)
    UU.allocate(1,lengthU);
    for (int i = 1; i <= lengthU; i++)
       UU(i) = i;
    ss = 1.0/dt;
    TTRACE(dt,ss);
    obs_catch.allocate(1,ngear,1,ntime);
    residuals.allocate(1,ntime,1,2*ngear+3);
    residuals.initialize();
    TRACE(catch_dat_file);
    cifstream cc(catch_dat_file);
    if (!cc)
    {
       clogf << "error opening " << catch_dat_file << endl;
       cerr  << "error opening " << catch_dat_file << endl;
       ad_exit(1);
    }

    clogf << "Opened " << catch_dat_file << endl;
    cc >> obs_catch;
    if (!cc)
    {
       clogf << "error reading " << catch_dat_file << endl;
       cerr  << "error reading " << catch_dat_file << endl;
       ad_exit(1);
    }
    ZeroCatch = 1.0; //1.0e-8;
    logZeroCatch = log(ZeroCatch);
    TTRACE(ZeroCatch,logZeroCatch)
    //ZeroF = 1.0e-11;
    //logZeroF = log(ZeroF);
    //TTRACE(ZeroF,logZeroF)
    //for (int g = 1; g <= ngear; g++)
    //  for (int t = 1; t <= ntime; t++)
    //     if (obs_catch(g,t) <= 0.0)
    //        obs_catch(g,t) = ZeroCatch;
    obs_catch = log(obs_catch+ZeroCatch);

    TRACE(forcing_dat_file)
    forcing_matrix.allocate(1,9,1,ntime);
    cifstream ff(forcing_dat_file);
    if (!ff)
    {
       clogf << "error opening " << forcing_dat_file << endl;
       cerr  << "error opening " << forcing_dat_file << endl;
       ad_exit(1);
    }
    clogf << "Opened " << forcing_dat_file << endl;
    ff >> forcing_matrix;
    if (!ff)
    {
       clogf << "error reading " << forcing_dat_file << endl;
       cerr  << "error reading " << forcing_dat_file << endl;
       ad_exit(1);
    }
    TRACE(fr);
    mean_immigrant_biomass = mean(forcing_matrix(fr));
    immigrant_biomass.allocate(1,ntime);
    //immigrant_biomass = forcing_matrix(fr);
    immigrant_biomass = mean_immigrant_biomass;
    TRACE(mean_immigrant_biomass)
    TRACE(immigrant_biomass)

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
  init_number logT12(phase_logT12);
  init_number logT21(phase_logT21);
  init_number logr(phase_logr);
  init_number logK(phase_logK);

  //init_vector logsdlogF(1,ngear,phase_logsdlogF);
  init_number logsdlogF(phase_logsdlogF);
  init_vector logsdlogPop(1,2,phase_logsdlogPop);
  //init_number logsdlogPop(phase_logsdlogPop);
  init_vector logsdlogYield(1,ngear,phase_logsdlogYield);

  // logit transformed porportion local
  init_number LmeanProportion_local(phase_LmeanProportion_local);
  init_number logsdLProportion_local(phase_logsdLProportion_local);

  // robust F likelihood
  init_vector Lpfat(1,ngear,phase_pfat);

  random_effects_vector U(1,lengthU);

  objective_function_value nll;

PRELIMINARY_CALCS_SECTION
    userfun_entries = 0;
    pininit = fexists(adstring(argv[0])+".pin");
    TRACE(pininit)
    if (!pininit)
    {
       logT12 = init_logT12;
       logT21 = init_logT21;
       TTRACE(logT12,logT21)
       logr = init_logr;
       TTRACE(logr,mfexp(logr))
       logK = init_logK;
       TTRACE(logK,mfexp(logK))

       logsdlogF = init_logsdlogF;
       for (int g = 1; g <= ngear; g++)
       {
          logsdlogYield(g) = init_logsdlogYield(g);
       }
       logsdlogPop = init_logsdlogPop;


       LmeanProportion_local = init_LmeanProportion_local;
       logsdLProportion_local = init_logsdLProportion_local;
       TTRACE(LmeanProportion_local,logsdLProportion_local)
       //double prop = 1.0/(1.0+mfexp(-value(LmeanProportion_local)));
       double prop = alogit(value(LmeanProportion_local));
       TTRACE(LmeanProportion_local,prop)

       if (!use_robustF)
       {
          phase_pfat = -1;
          init_pfat = 1e-25;
       }
       for (int g = 1; g <= ngear; g++)
       {
          Lpfat(g) = logit(init_pfat(g));
       }
       TRACE(init_pfat)
       TRACE(Lpfat)

       double K = mfexp(value(logK));
       double Pop1 = prop*K;
       double Pop2 = K-Pop1;
       TTRACE(Pop1,Pop2)
       TTRACE(log(Pop1),log(Pop2))

       dmatrix Ferr(1,ntime,1,ngear); Ferr.fill_randn(77);
       dmatrix logPop1Err(1,ntime,1,2); logPop1Err.fill_randn(79);
       dmatrix logPop2Err(1,ntime,1,2); logPop2Err.fill_randn(75);
       int ut = 0;
       TRACE(ut)
       for (int t = 1; t <= ntime; t++)
       {   
          for (int g = 1; g <= ngear; g++)
          {
          // U(++ut) = log(Fmsy)+0.1*Ferr(t,g); // Fmsy
          // U(++ut) = log(1.0e-8+1.0e-7*Ferr(t,g));
             U(++ut) = -18*mfexp(0.001*Ferr(t,g));
             if ((t<=2) && (g==1))
               TTRACE(Ferr(t,g),U(ut))
          }
       }
       TTRACE(ut,utPop1)
       for (int t = 1; t <= ntime; t++)
       {   
          //logPop(t,1) = log(Pop1)+0.5*logPop1Err(t,1);
          U(++ut) = log(Pop1)+0.5*logPop1Err(t,1);
       }
       TTRACE(ut,utPop2)
       for (int t = 1; t <= ntime; t++)
       {   
          //logPop(t,2) = log(Pop2)+0.5*logPop2Err(t,2);
          U(++ut) = log(Pop2)+0.5*logPop2Err(t,2);
       }
       TRACE(ut)
       ofstream pin("xssams.p00");
       if (!pin)
       {
          cerr << "Error creating xssams.p00" << endl;
          ad_exit(1);
       }
       PINOUT(logT12)
       PINOUT(logT21)
       PINOUT(logr)
       PINOUT(logK)
       PINOUT(logsdlogF)
       PINOUT(logsdlogPop)
       PINOUT(logsdlogYield)
       PINOUT(LmeanProportion_local)
       PINOUT(logsdLProportion_local)
       PINOUT(alogit(value(Lpfat)))
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
          cerr << "Error writing xssams.p00" << endl;
          ad_exit(1);
       }
       else
       {
          clogf << "Successfully created xssams.p00" << endl;
          cout << "Successfully created xssams.p00" << endl;
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
    TRACE(logr)
    TRACE(logK)
    TRACE(logsdlogF)
    TRACE(logsdlogPop)
    TRACE(logsdlogYield)
    TRACE(LmeanProportion_local)
    TRACE(logsdLProportion_local)
    TRACE(alogit(value(Lpfat)))
    TRACE(U)
    TTRACE(U(utPop1+1),U(utPop2+1))
    TRACE(obs_catch)
    TRACE (trace_init_pars)

    clogf << endl;
    //if (1) ad_exit(1);
  

RUNTIME_SECTION
  convergence_criteria .1, 1e-5


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
    TRACE(logr)
    TRACE(logK)
    TRACE(logsdlogF)
    TRACE(logsdlogPop)
    TRACE(logsdlogYield)
    TRACE(LmeanProportion_local)
    TRACE(logsdLProportion_local)
    TRACE(U)
    TTRACE(U(utPop1+1),U(utPop2+1))
    TRACE(obs_catch)
    trace_init_pars = 0;
    TRACE (trace_init_pars)
    clogf << endl;
    //if (1) ad_exit(1);
  }

  nll = 0.0;

  for (int t = 2; t <= maxtime; t++)
  {
     step(t, U(Fndxl(t-1),Fndxu(t-1)), U(Fndxl(t),Fndxu(t)), logsdlogF,
             U(utPop1+t-1), U(utPop1+t), U(utPop2+t-1), U(utPop2+t), logsdlogPop,
             logr,logK,logT12,logT21,LmeanProportion_local,logsdLProportion_local,Lpfat);

  }

  for (int t = 1; t <= maxtime; t++)
  {
     obs(t,U(Fndxl(t),Fndxu(t)),U(utPop1+t-1),U(utPop1+t),U(utPop2+t-1),U(utPop2+t),logsdlogYield);
  }

  //TRACE(++userfun_entries)
  ++userfun_entries;
  int status_print = ntime;
  if (userfun_entries > lengthU)
     status_print = lengthU;
  if ((userfun_entries == 1) || (userfun_entries % status_print) == 0)
  {
     write_status(clogf);
  }
  /*
  clogf << "###param";
  clogf << " " << userfun_entries;
  clogf << " " << nll;
  clogf << " " << logT12;
  clogf << " " << logT21;
  clogf << " " << logr;
  clogf << " " << logK;
  clogf << " " << logsdlogF;
  clogf << " " << logsdlogPop;
  clogf << " " << logsdlogYield;
  clogf << " " << LmeanProportion_local;
  clogf << " " << logsdLProportion_local;
  clogf << " " << Lpfat;
  clogf << endl;
  */

  //if (1) ad_exit(1);


SEPARABLE_FUNCTION void step(const int t, const dvar_vector& f1, const dvar_vector& f2, const dvariable& lsdlogF, const dvariable& p11, const dvariable p12, const dvariable& p21, const dvariable p22, const dvar_vector& lsdlogPop, const dvariable& lr, const dvariable& lK, const dvariable& lT12, const dvariable& lT21, const dvariable& LmPropL, const dvariable& lsdLProportion_local, const dvar_vector& Lpf)
  // p11 U(utPop1+t-1) log N1 at start of time step
  // p12 U(utPop1+t)   log N1 at end   of time step
  // p21 U(utPop2+t-1) log N2 at start of time step
  // p22 U(utPop2+t)   log N2 at end   of time step

  dvariable varlogF = square(mfexp(lsdlogF));
  dvar_vector varlogPop = square(mfexp(lsdlogPop));
  dvariable r = mfexp(lr);
  dvariable K = mfexp(lK);
  dvariable T12 =mfexp(lT12);
  dvariable T21 =mfexp(lT21);
  dvariable LmeanPropL = LmPropL;
  dvariable varLPropL = square(mfexp(lsdLProportion_local));
  dvariable nextLogN1;
  dvariable nextLogN2;
  dvar_vector ft1(1,ngear);
  dvar_vector ft2(1,ngear);
  for (int g = 1; g <= ngear; g++)
  {
     ft1(g) = f1(f1.indexmin()+g-1);
     ft2(g) = f2(f2.indexmin()+g-1);
  }

  dvariable Fnll = 0.0;
  //if (use_robustF)
  //   dvar_vector pfat(1,ngear);
  for (int g = 1; g <= ngear; g++)
  {
     if (use_robustF == 1)  // normal + t distribution
     {
        dvariable z = square(ft1(g)-ft2(g))/varlogF;
        dvariable norm_part = sqrt(varlogF) + 0.5*LOG_TWO_M_PI + z;
        dvariable fat_part = LOG_M_PI + log(1.0 + z);
        dvariable pfat = alogit(Lpf(g));
        Fnll += log((1.0-pfat)*mfexp(norm_part) + pfat*mfexp(fat_part));
        //TTRACE(norm_part,fat_part)
	//TTRACE(pfat,Lpf(g))
     }
     else if (use_robustF == 2) // from newreg2.cpp
     {
        double width=3.0;
        double width2=width*width;
        const double alpha = 0.7;
        const double a2 = square(alpha);
        dvariable diff2 = square(ft1(g)-ft2(g));
        dvariable v_hat = diff2+1.0e-80;

        dvariable pcon = alogit(Lpf(g));
        dvariable b=2.*pcon/(width*sqrt(PI));  // This is the weight for the "robustifying" term

        dvariable norm_part = log(1.0-pcon)*mfexp(-diff2/(2.0*a2*v_hat));
   
        dvariable fat_part =  b/(1.+pow(diff2/(width2*a2*v_hat),2));
        //TTRACE(norm_part,fat_part)
        //TRACE(mfexp(norm_part + fat_part))
        Fnll += norm_part + fat_part;
     }

     else
     {
        Fnll += 0.5*(log(TWO_M_PI*varlogF) + square(ft1(g)-ft2(g))/varlogF);
     }
     //TRACE(Fnll)
     if (isnan(value(Fnll)))
     {
        TTRACE(Fnll,varlogF)
        TTRACE(ft1(g),ft2(g))
        TTRACE(t,g)
        write_status(clogf);
        ad_exit(1);
     }
  } 

  //dvariable sumFg = 0.0; // total fishing mortality
  dvariable sumFg = sum(mfexp(ft1)); // total fishing mortality
  nextLogN1 = p11;
  nextLogN2 = p21;
  dvariable prevLogN1;
  dvariable prevLogN2;
  dvariable prevN1;
  dvariable prevN2;
  for (int s = 1; s <= ss; s++)
  {
     //TTRACE(s,dt)
     prevLogN1 = nextLogN1;
     prevLogN2 = nextLogN2;
     //TTRACE(prevLogN1,prevLogN2)
     prevN1 = mfexp(prevLogN1);
     prevN2 = mfexp(prevLogN2);
     //TTRACE(prevN1,prevN2)

     nextLogN1 += dt*(r*(1.0 - prevN1/K) - sumFg - T12 - r*prevN2/K);
     nextLogN2 += dt*(r*(1.0 - prevN2/K) - sumFg - T12 - r*prevN1/K + T21*immigrant_biomass(t)/prevN2);
     //TTRACE(nextLogN1,nextLogN2)
  }
  

  if ( isnan(value(nextLogN1)) || isnan(value(nextLogN2)) ||
       isinf(value(nextLogN1)) || isinf(value(nextLogN2)) )
  {
     TTRACE(nextLogN1,nextLogN2)
     TTRACE(prevLogN1,prevLogN2)
     TTRACE(prevN1,prevN2)
     TTRACE(r,K)
     write_status(clogf);
     ad_exit(1);
  }

  dvariable Pnll = 0.0;
  dvariable diff = log(mfexp(p12)+mfexp(p22))-log(mfexp(nextLogN1)+mfexp(nextLogN2));
  // process error N.
  //Pnll += 0.5*(log(TWO_M_PI*varlogPop) + square(diff)/varlogPop);
  //if (isnan(value(Pnll)))
  //{
  //   TRACE(Pnll)
  //   write_status(clogf);
  //   ad_exit(1);
  //}

  // process error N1
  Pnll += 0.5*(log(TWO_M_PI*varlogPop(1)) + square(p12-nextLogN1)/varlogPop(1));
  if (isnan(value(Pnll)))
  {
     TRACE(Pnll)
     write_status(clogf);
     ad_exit(1);
  }

  // process error N2
  Pnll += 0.5*(log(TWO_M_PI*varlogPop(2)) + square(p22-nextLogN2)/varlogPop(2));
  if (isnan(value(Pnll)))
  {
     TRACE(Pnll)
     write_status(clogf);
     ad_exit(1);
  }

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

  

SEPARABLE_FUNCTION void obs(const int t, const dvar_vector& f,const dvariable& pop11, const dvariable& pop21,  const dvariable& pop12, const dvariable& pop22, const dvar_vector& logsdlogYield)
  // pop11 U(utPop1+t-1)
  // pop21 U(utPop1+t)
  // pop12 U(utPop2+t-1)
  // pop22 U(utPop2+t)

  dvar_vector ft(1,ngear); ///< log fishing mortality by gear
  for (int g = 1; g <= ngear; g++)
  {
     //TTRACE(g,UU(f.indexmin()+g-1))
     ft(g) = f(f.indexmin()+g-1);
  }

  dvar_vector log_pred_yield(1,ngear);
  // sum of the average populations sizes over time step
  dvariable log_total_mean_pop = log( 0.5*(mfexp(pop11) + mfexp(pop21) +   // population 1
                                           mfexp(pop12) + mfexp(pop22)) ); // population 2
  for (int g = 1; g <= ngear; g++)
  {
     log_pred_yield(g) =  ft(g) + log_total_mean_pop;
  }

  dvar_vector varlogYield = square(mfexp(logsdlogYield));
  dvariable Ynll = 0.0;
  for (int g = 1; g <= ngear; g++)
  {
        // observation error
        Ynll += 0.5*(log(TWO_M_PI*varlogYield(g)) + square(obs_catch(g,t)-log_pred_yield(g))/varlogYield(g));
        if (isnan(value(Ynll)))
        {
           TRACE(Ynll)
           TTRACE(t,g)
           write_status(clogf);
           ad_exit(1);
        }
  }
  //clogf << t << " " << Ynll << " " << nll << endl;
  nll += Ynll;

  residuals(t,1) = value(pop21);
  residuals(t,2) = value(pop22);
  double propLa = value(pop21) - value(pop22);
  //double propLb = mfexp(value(pop21))/(mfexp(value(pop21)) + mfexp(value(pop22)));
  //TTRACE(alogit(propLa),propLb)
  residuals(t,3) = alogit(propLa);

  int rc = 3 ;
  for (int g = 1; g <= ngear; g++)
    residuals(t, ++rc) = value(ft(g));
  for (int g = 1; g <= ngear; g++)
    residuals(t, ++rc) = value(log_pred_yield(g));

FUNCTION void write_status(ofstream& s)
    //double prop = 1.0/(1.0+mfexp(value(-LmeanProportion_local)));
    double prop = alogit(value(LmeanProportion_local));
    s << "\n# Status after "<< userfun_entries << " PROCEDURE_SECTION entries;" << endl;
    s << "# nll = " << value(nll) << endl;
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
    s << "# pfat = " << alogit(value(Lpfat)) << " (" << active(Lpfat) <<")" << endl;
    s << "# Residuals:" << endl;
    s << "  t    pop1   pop2  propL";
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
          s << " " << obs_catch(g,t);
       s << endl;
    }



REPORT_SECTION
    //double prop = 1.0/(1.0+mfexp(value(-LmeanProportion_local)));
    double prop = alogit(value(-LmeanProportion_local));
    REPORT(logT21)
    REPORT(logT12)
    REPORT(logr)
    REPORT(logK)
    REPORT(LmeanProportion_local)
    REPORT(logsdLProportion_local)
    REPORT(prop)
    REPORT(logsdlogF)
    REPORT(logsdlogPop)
    REPORT(logsdlogYield)
    write_status(report);

