GLOBALS_SECTION;
  #include <math.h>
  #include <adstring.hpp>
  #include "trace.h"
  //#include <df1b2fun.h>

  //#include "nLogNormal.h"
 
  #undef PINOUT
  #define PINOUT(object) pin << "# " << #object ":\n" << setprecision(5) << object << endl;
  #undef REPORT
  #define REPORT(object) report << "# " << #object ":\n" << setprecision(5) << object << endl;

  ofstream clogf;
  const double TWO_M_PI = 2.0*M_PI;

  // ./xssams -noinit -nr 2 -est -l2 10000000  -l3 10000000
  int fexists(const adstring& filename)
  {
    std::ifstream ifile(filename);
    if (!ifile)
       return 0;
    else
       return 1;
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

  //pad();

DATA_SECTION;
  init_adstring catch_dat_file;
  init_adstring forcing_dat_file;

  init_int ntime;
  init_int ngear;
  init_number dt; // integration time step
  init_int fr;
  number ss;      // 1/dt number of interations in the integration
  //init_adstring catch_dat_file(pxml);
  matrix obs_catch;
  //matrix log_obs_catch;
  matrix forcing_matrix;
  vector immigrant_biomass;
  int pininit;
  number Fmsy;
  number Bmsy;
  number MSY;
  int lengthU;
  int utPop1;
  int utPop2;
  int userfun_entries;
  int maxtime;
  ivector Fndxl;
  ivector Fndxu;
  int trace_init_pars;
  matrix residuals;

 LOCAL_CALCS;
    maxtime = ntime;
    TTRACE(ntime,maxtime);
    TRACE(ngear);
    lengthU = ntime*(ngear+2);
    TRACE(lengthU)
    ss = 1.0/dt;
    TTRACE(dt,ss);
    obs_catch.allocate(1,ngear,1,ntime);
    residuals.allocate(1,ntime,1,2*ngear+2);
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
    for (int g = 1; g <= ngear; g++)
      for (int t = 1; t <= ntime; t++)
         if (obs_catch(g,t) <= 0.0)
            obs_catch(g,t) = 1.0e-8;
    obs_catch = log(obs_catch);

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
    immigrant_biomass = forcing_matrix(fr);
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

PARAMETER_SECTION;
  init_number logT12(-1);
  init_number logT21(-1);
  number T12;
  number T21;
  init_bounded_number r(0.0,0.25,1);
  init_bounded_number K(1.0e-8,1e8,-1);

  init_bounded_vector logsdlogF(1,ngear,-5,5,2);
  //init_vector logsdlogF(1,ngear,1);

  init_bounded_vector logsdlogPop(1,2,-5,5,2);
  //init_vector logsdlogPop(1,2,1);

  init_bounded_vector logsdlogYield(1,ngear,-5,5,1);
  //init_vector logsdlogYield(1,ngear,1);

  init_number LmeanProportion_local(2); // logit transformed porportion local
  init_number logsdLProportion_local(1);

  random_effects_vector U(1,lengthU,1);

  objective_function_value nll;

PRELIMINARY_CALCS_SECTION
    userfun_entries = 0;
    pininit = fexists(adstring(argv[0])+".pin");
    TRACE(pininit)
    if (!pininit)
    {
       logT12 = -8.0;
       logT21 = -8.0;;
       TTRACE(logT12,logT21)
       r = 0.15;
       //Fmsy = 0.5*value(r);
       //TTRACE(r,Fmsy)
       //MSY = 5000.0;
       //Bmsy = MSY/Fmsy;
       //K = r*Bmsy/(r-Fmsy);
       //TTRACE(Bmsy,K)
       K = 1e6;

       for (int g = 1; g <= ngear; g++)
       {
          logsdlogF(g) = 1.1;
          logsdlogYield(g) = 1.1;
       }
       logsdlogPop(1) = 1.1;
       logsdlogPop(2) = 1.1;

       double prop = 0.9;
       LmeanProportion_local = log(prop/(1.0-prop));
       logsdLProportion_local = 1.5;

       double Pop1 = value(prop*K);
       double Pop2 = value(K-Pop1);
       TTRACE(Pop1,Pop2)
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
             U(++ut) = -18*exp(0.001*Ferr(t,g));
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
       PINOUT(r)
       PINOUT(K)
       PINOUT(logsdlogF)
       PINOUT(logsdlogPop)
       PINOUT(logsdlogYield)
       PINOUT(LmeanProportion_local)
       PINOUT(logsdLProportion_local)
       //PINOUT(U)
       pin << "# U:" << endl;
       pin << "#   F(t,g):" << endl;
       for (int tg = 1; tg <= (ntime*ngear); tg++)
           pin << " " << U(tg);
       pin << "\n#   Pop1(t):"<< endl;
       //int tt1 = 1 + ntime*ngear;
       //int tt2 = tt1 + ntime;
       //TTRACE(tt1,tt2)
       //for (int tt = tt1; tt <= tt2; tt++)
       for (int tt = 1; tt <= ntime; tt++)
          pin << " " << U(utPop1+tt);
       pin << "\n#   Pop2(t):"<< endl;
       //tt1 = tt2 + 1;
       //tt2 = tt1 + ntime;
       //TTRACE(tt1,tt2)
       //for (int tt = tt1; tt <= tt2; tt++)
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
          cout << "Successfully created xssams.p00" << endl;
          clogf << "Successfully created xssams.p00" << endl;
       }
    }
    trace_init_pars = 1;

    clogf << "\nAt end of PRELIMINARY_CALCS_SECTION:"<<endl;
    TRACE(logT12)
    TRACE(logT21)
    TRACE(r)
    TRACE(K)
    TRACE(logsdlogF)
    TRACE(logsdlogPop)
    TRACE(logsdlogYield)
    TRACE(U)
    TTRACE(U(utPop1+1),U(utPop2+1))
    TRACE(obs_catch)
    clogf << endl;
    //if (1) ad_exit(1);


PROCEDURE_SECTION
  if (trace_init_pars)
  {
    clogf << "\nInitial step in to PROCEDURE_SECTION:" << endl;
    TRACE(logT12)
    TRACE(logT21)
    TRACE(r)
    TRACE(K)
    TRACE(logsdlogF)
    TRACE(logsdlogPop)
    TRACE(logsdlogYield)
    TRACE(U)
    TRACE(obs_catch)
    trace_init_pars = 0;
    TRACE (trace_init_pars)
    clogf << endl;
    //if (1) ad_exit(1);
  }

  nll = 0.0;
  T12 = exp(logT12);
  T21 = exp(logT21);

  //for (int t = 2; t <= ntime; t++)
  for (int t = 2; t <= maxtime; t++)
  {
     step(t,U(Fndxl(t-1),Fndxu(t-1)),U(Fndxl(t),Fndxu(t)),logsdlogF,U(utPop1+t-1),U(utPop1+t),U(utPop2+t-1),U(utPop2+t),logsdlogPop,r,K,T12,T21,LmeanProportion_local,logsdLProportion_local);
  }

  //for (int t = 1; t <= ntime; t++)
  for (int t = 1; t <= maxtime; t++)
  {
     //TRACE(t)
     obs(t,U(Fndxl(t),Fndxu(t)),U(utPop1+t),U(utPop2+t),logsdlogYield);
  }

  TRACE(++userfun_entries)
  if ((userfun_entries % ntime) == 0)
  {
     write_status(clogf);
  }


SEPARABLE_FUNCTION void step(const int t, const dvar_vector& f1, const dvar_vector& f2, const dvar_vector& logsdF, const dvariable& p11, const dvariable p12, const dvariable& p21, const dvariable p22, const dvar_vector& logsdPop, const dvariable& r, const dvariable& K, const dvariable& T12, const dvariable& T21, const dvariable& LmeanPropL, const dvariable& logsdLPropL)

  dvar_vector ft1(1,ngear);
  dvar_vector ft2(1,ngear);
  for (int g = 1; g <= ngear; g++)
  {
     ft1(g) = f1(f1.indexmin()+g-1);
     ft2(g) = f2(f2.indexmin()+g-1);
  }
  if (t == 2)
  {
     dvariable logF0 = log(1.0e-8);
     for (int g = 1; g <= ngear; g++)
     {
        dvariable varlogF = square(exp(logsdF(g)));
        nll += 0.5*(log(TWO_M_PI*varlogF) + square(ft1(g)-logF0)/varlogF);
     }
  }

  dvariable sumFg = sum(ft1); // total fishing mortality

  for (int g = 1; g <= ngear; g++)
  {
     dvariable varlogF = square(exp(logsdF(g)));
     // fishing mortality prior
     nll += 0.5*(log(TWO_M_PI*varlogF) + square(ft1(g)-ft2(g))/varlogF);
  } 

  dvar_vector varPop = square(exp(logsdPop));
  dvariable nextLogN1;
  dvariable nextLogN2;
  if (t == 2)
  {
     dvariable tprop = 1.0/(1.0+exp(-LmeanPropL));
     //TRACE(tprop)
     nextLogN1 = log(tprop*K);
     nextLogN2 = log(K - tprop*K);
     //TTRACE(nextLogN1,nextLogN2)
     //TTRACE(exp(nextLogN1),exp(nextLogN2))
     nll += 0.5*(log(TWO_M_PI*varPop(1)) + square(nextLogN1-p11)/varPop(1));
     nll += 0.5*(log(TWO_M_PI*varPop(2)) + square(nextLogN2-p21)/varPop(2));
  }
  
  nextLogN1 = p11;
  nextLogN2 = p21;
  //TTRACE(exp(nextLogN1),exp(nextLogN2))
  dvariable prevLogN1;
  dvariable prevLogN2;
  dvariable prevN1;
  dvariable prevN2;
  for (int s = 1; s <= ss; s++)
  {
     prevLogN1 = nextLogN1;
     prevLogN2 = nextLogN2;
     prevN1 = exp(prevLogN1);
     prevN2 = exp(prevLogN2);

     nextLogN1 += dt*(r*(1.0 - prevN1/K) - sumFg - T12 - r*prevN2/K);
     nextLogN2 += dt*(r*(1.0 - prevN2/K) - sumFg - T12 - r*prevN1/K + T21*immigrant_biomass(t)/prevN2);
  }
  // process error
  nll += 0.5*(log(TWO_M_PI*varPop(1)) + square(nextLogN1-p12)/varPop(1));
  nll += 0.5*(log(TWO_M_PI*varPop(2)) + square(nextLogN2-p22)/varPop(2));

  // proportion local
  // Proportion_local = exp(nextLogN1)/(exp(nextLogN1)+exp(nextLogN2));
  dvariable LpropL = nextLogN1 - nextLogN2;
  dvariable varLPropL = square(exp(logsdLPropL));
  nll += 0.5*(log(TWO_M_PI*varLPropL) + square(LpropL - LmeanPropL)/varLPropL);

  

SEPARABLE_FUNCTION void obs(const int t, const dvar_vector& f, const dvariable& pop1, const dvariable& pop2, const dvar_vector logsdYield)
  dvar_vector ft(1,ngear);
  for (int g = 1; g <= ngear; g++)
  {
     ft(g) = f(f.indexmin()+g-1);
  }

  dvar_vector pred_yield(1,ngear);
  dvariable logTotalPop = pop1 + pop2;
  for (int g = 1; g <= ngear; g++)
  {
     pred_yield(g) =  ft(g) + logTotalPop;
  }

  dvariable Ynll = 0.0;
  for (int g = 1; g <= ngear; g++)
  {
     if (obs_catch(g,t) > -18.0)
     {
     dvariable varYield = square(exp(logsdYield(g)));
     // observation error
     Ynll += 0.5*(log(TWO_M_PI*varYield) + square(obs_catch(g,t)-pred_yield(g))/varYield);
     }
  }
  nll += Ynll;

  residuals(t,1) = value(pop1);
  residuals(t,2) = value(pop2);
  int rc = 2 ;
  for (int g = 1; g <= ngear; g++)
    residuals(t, ++rc) = value(ft(g));
  for (int g = 1; g <= ngear; g++)
    residuals(t, ++rc) = value(pred_yield(g));

FUNCTION void write_status(ofstream& s)
    double prop = 1.0/(1.0+exp(value(-LmeanProportion_local)));
    s << "\n# Status after "<< userfun_entries << " PROCEDURE_SECTION entries;" << endl;
    s << "# nll = " << value(nll) << endl;
    s << "# nvar = " << initial_params::nvarcalc() << endl;
    s << "#  T12 = " << T12 << endl;
    s << "#  T21 = " << T21 << endl;
    s << "#    r = " << r << endl;
    s << "#    K = " << K << endl;
    s << "# LmeanProportion_local = " << LmeanProportion_local << endl;
    s << "#                  prop = " << prop << endl;
    s << "# logsdLProportion_local = " << logsdLProportion_local << endl;
    s << "#     logsdlogF: " << logsdlogF << endl;
    s << "#   logsdlogPop: " << logsdlogPop << endl;
    s << "# logsdlogYield: " << logsdlogYield << endl;
    s << "# Residuals:" << endl;
    s << "  t    pop1   pop2";
    for (int g = 1; g <= ngear; g++)
       s << "     F" << g;
    for (int g = 1; g <= ngear; g++)
       s << "  predC" << g;
    for (int g = 1; g <= ngear; g++)
       s << "   obsC" << g;
    s << endl;
    for (int t = 1; t <= ntime; t++)
    {
       s << t << " " << residuals(t,1) << " " << residuals(t,2);
       int rc = 2;
       for (int g = 1; g <= 2*ngear; g++)
          s << " " << residuals(t,++rc);
       for (int g = 1; g <= ngear; g++)
          s << " " << obs_catch(g,t);
       s << endl;
    }



REPORT_SECTION
    double prop = 1.0/(1.0+exp(value(-LmeanProportion_local)));
    REPORT(T21)
    REPORT(T12)
    REPORT(r)
    REPORT(K)
    REPORT(LmeanProportion_local)
    REPORT(logsdLProportion_local)
    REPORT(prop)
    REPORT(logsdlogF)
    REPORT(logsdlogPop)
    REPORT(logsdlogYield)
    write_status(report);

