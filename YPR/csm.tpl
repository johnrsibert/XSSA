// Modification of csm.tpl originally written in 2001 by Shiham Adam.
// It is a simple tag attrition model based on fadmodel9.tpl intended 
// to estimate fishing and natural mortality for use in yield per recruit analysi

GLOBALS_SECTION
  #include <sstream>
  #include <fvar.hpp>
  #include "trace.h"
 
  #undef REPORT
  #define REPORT(object) report << "#" << #object"\n" << setprecision(5) << object << endl;

  ofstream clogf;
  int ndx = 1;
  int sze = 2;
  int num = 3;

  void strat(const int max_days, const int bin_days) 
  {
     using namespace std;
     stringstream cs;

     cs << "strat " << max_days << " " << bin_days << ends; 
     int rv = system(cs.str().c_str());
     if (rv)
     {
       cerr << "Error attempting to spawn '" << cs.str() << "'" << endl;
       exit(1);
     }
     cout << "\nSpawned '" << cs.str() << "'" << endl; 
   }
  
TOP_OF_MAIN_SECTION
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(4000000);
  arrmblsize = 3000000;

  adstring logname(adstring(argv[0])+"_program.log");
  clogf.open(logname);
  if ( !clogf ) {
    cerr << "Cannot open program log " << logname << endl;
    ad_exit(1);
  }
  cout << "Opened program log: " << logname << endl;

DATA_SECTION
  init_int max_days;
  !! TRACE(max_days);
  init_int bin_days;  
  !! TRACE(bin_days);
  init_int min_nos_in_cohort;
  !! TRACE(min_nos_in_cohort);
  init_int nexclude; //no of steps to eliminate
  !! TRACE(nexclude); 
  init_number T1shed;
  !! TRACE(T1shed);
  //init_int nstep;      
  //!! TRACE(nstep); 
  init_int ngroup_MF;  
  !! TRACE(ngroup_MF);                       
  init_ivector MF_sz_grp(1,ngroup_MF);  
  !! TRACE(MF_sz_grp);     
  init_number M_curve_pen_weight;
  !! TRACE(M_curve_pen_weight);

  init_int st;
  !! TRACE(st); 
  init_int fn;
  !! TRACE(fn);
   
  !! strat(max_days, bin_days);      //send the conosle command
  !! ad_comm::change_datafile_name("yft.dat");     //swap file
  //!! ad_comm::change_datafile_name("bench.dat");     //swap file
  init_int min_size;
  init_int max_size;
  !! TTRACE(min_size, max_size);
  init_int ncohort; 
  !! TRACE(ncohort);
  init_imatrix rel_data(1,ncohort,1,3);
  !! TRACE(rel_data(ncohort))
  init_int ntime;
  !! TRACE(ntime);
  init_imatrix obs_caps(1,ncohort,0,ntime);
  !! TRACE(obs_caps(ncohort))
  ivector fish_size(min_size,2*max_size);               //set up size vector  
  number dt;
  //number tt;
  vector nobs(1,ngroup_MF);

  
PARAMETER_SECTION
  init_bounded_vector F(1,ngroup_MF,0.0,0.4);
  init_bounded_vector M(1,ngroup_MF,0.0,0.4,2);

  sdreport_vector qF(1,ngroup_MF);
  sdreport_vector qM(1,ngroup_MF);

  number Z;
  number N1;
  number N2;
  number numrel;
  vector Fdum(0,nexclude); 
  number prop_remain;   
  matrix pred_caps(1,ncohort,0,ntime);
  number M_curve_penalty;
  
  number ll;
  objective_function_value L;

PRELIMINARY_CALCS_SECTION
  //sets up the size index vector 
  fish_size = ngroup_MF; // initialize to index of largest group
  int ss = min_size;
  for(int a = 1; (a < ngroup_MF) && (ss <= max_size); a++)
  {
    while( ss <= MF_sz_grp(a+1))
    {
       fish_size(ss) = a;
       ss++;
    }
  }
  TRACE(fish_size); 

  //calculate dt
  //dt = double(bin_days)/double(nstep);
  //dt = double(1)/double(nstep);
  dt = double(bin_days);
  TRACE(dt);

  

RUNTIME_SECTION
  convergence_criteria 1.e-2, 1.e-3, 1.e-4;
  maximum_function_evaluations 1500, 1500, 800;

PROCEDURE_SECTION
  pred_caps.initialize();
  L = 0.0;
  ll = 0.0;
  nobs.initialize();

  //set up curvature penalty on M
  dvariable pen_sum = 0.0;
  for(int k = 2; k < ngroup_MF; k++ ) 
  {
     pen_sum += M(k-1) - 2.0 * M(k) + M(k+1);
  }
  M_curve_penalty = pen_sum * M_curve_pen_weight;

  //int c = 156;
  for (int c = st; c <= fn;  c++)
  {
     //N1 = rel_data(c,num);  // nos released in the cohort 
     //tag shedding will be accounted for cohorts > min_nos_in_cohort! 
     N1 = (rel_data(c,num) > min_nos_in_cohort) ? T1shed * rel_data(c,num) :  rel_data(c,num);

     double rel_len = rel_data(c,sze);  //get release length
     double midtime; 
     double pred_len;
     int ipl; 
     int ndx;

     for (int t = 0; t <= ntime; t++)  
     {
        midtime = bin_days * (double(t) + 0.5); //TRACE(midtime); 
        // BET
        //pred_len =  (181.7 - rel_len) * (1 - mfexp(-0.00068767* midtime)) + rel_len;      
        //YFT
        pred_len =  (166.4 - rel_len) * (1 - mfexp(-0.0006887* midtime)) + rel_len;      
        //TTRACE(rel_len,pred_len)
        ipl = (pred_len + 0.5);
        ndx = fish_size(ipl);
        nobs(ndx) ++;

        Z = F(ndx)+M(ndx);
        N2 = N1 * mfexp( (-Z*dt)   );
        if (isnan(value(N2)))
        {
           TTRACE(N1,N2)
           TTRACE(c,t)
           ad_exit(1);
        }
        pred_caps(c, t) = N1 * F(ndx)/Z *(1.0-mfexp(-Z*dt)  );
        if (isnan(value(pred_caps(c,t))))
        {
           TTRACE(N1,N2)
           TTRACE(c,t)
           ad_exit(1);
        }
        N1 = N2;


     } //for int t

   } //for cohort
  
  //Poission Likelihood function
  for(int i = st; i <= fn; i++)
  {
     for(int j = 0; j <= ntime; j++)
     { 
        //L -= obs_caps(i,j) * log(sfabs(pred_caps(i,j))+1e-20) - pred_caps(i,j);
        dvariable testL = obs_caps(i,j) * log(sfabs(pred_caps(i,j))+1e-20) - pred_caps(i,j);
        if (isnan(value(testL)))
        {
           TRACE(testL)
           TTRACE(i,j)
           TTRACE(obs_caps(i,j),pred_caps(i,j))
           TRACE(F)
           TRACE(M)
        }
        else
           L -= testL;
     }
  } 
  //TRACE (L);
  
  if(active(M))
  {
     //TTRACE(L,M_curve_penalty)
     L -= M_curve_penalty;
  } 
  qM = 91.25*M;
  qF = 91.25*F;
  


FUNCTION double CalcF(double rec, int rel, double tM) 
  { 
   /********************************************************************
   * this function iteratively calculates the value of the special F   *
   * that would make the observed and predicted reacptures the same    *
   ********************************************************************/
  //TRACE(rec);
  //TRACE(rel);

  double Z;
  double EZ;
  double RecHat;  //predicted recpature
  double dRecHat; //derivative of predicted recapture
  double delta;
  double dumF; 
  
  if(rec == 0.0)
  {
   dumF = 0.0;
   return (dumF);
   //exit(1);
  } 
  
  dumF = rec/rel;   //starting value 
  //dumF = 0.0001;   //starting value 
  for(int its = 1 ; its <= 100; its++)
  { 
     Z = dumF + tM;
     EZ = mfexp(-Z);
     RecHat = dumF / Z * (1 - EZ) * rel - rec;  
     dRecHat = rel * (Z * (1 - EZ + dumF * EZ) - dumF * (1 - EZ)) / Z / Z;
     delta = RecHat / dRecHat;
     dumF = dumF - delta;
  
     if ( fabs(delta) < 0.000001 && fabs(RecHat < 0.000001) )
      break; 
   }
   
  return(dumF);
  }


REPORT_SECTION
  REPORT(L)
  REPORT(M_curve_penalty)
  REPORT(F)
  REPORT(qF)
  REPORT(M)
  REPORT(qM)
  report << endl;
  report << "days" << "  " << "obs" << "   " << "pred" << endl; 
  for (int t = 0; t <= ntime; t++)
  {
    double obs = 0.0;
    double pre = 0.0;
    
    for (int cohort = st; cohort <= fn; cohort++)
     {
       obs += obs_caps(cohort,t);
       pre += value(pred_caps(cohort,t));
     }
     report <<  t* bin_days << "  "
           <<  obs         << "  "
           <<  pre         << endl;
  }
  report << endl;
  const double wa = 2.512e-05;
  const double wb = 2.9396;
  REPORT(ngroup_MF)
  report << "n len wt nobs M F" << endl;
  for (int n = 1; n <= ngroup_MF; n++)
  {
     double len = (double)MF_sz_grp(n) + 5.0;
     double wt = wa*pow(len,wb);
     report << n << " " << len << " " << wt << " " << nobs(n) << " " 
            << qM(n) << " " << qF(n) << endl;
  }
