#include <fvar.hpp>
#include <cifstrem.h>
#include "trace.h"
extern "C" {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}

int main(int argc, char* argv[])
{
 ad_exit=&ad_boundf;

 ofstream clogf("strat.log");
 TTRACE(argc, argv[0]);
 int narg = 2;
 if (argc <= narg)
 {
  cerr <<"Usage: strat  'max_day'  'bin_days' " << endl;
  exit(1); 
 }
 int max_days = atoi(argv[1]);
 int bin_days = atoi(argv[2]);
 TTRACE(max_days,bin_days);
 
 cifstream rel("yft_rel.txt");
 int ncohort;  rel >> ncohort;   TRACE(ncohort);
 imatrix rel_data(1,ncohort,1,3);
 rel >> rel_data;  //read data
 rel.close();
 clogf << rel_data(ncohort) << endl;
 
 cifstream rec("yft_rec.txt");
 int nrecap;  rec >> nrecap; TRACE(nrecap); 
 imatrix rec_data(1,nrecap,1,2);
 rec >> rec_data;  //read data
 rec.close();
 clogf << rec_data(nrecap) << endl;
 

 //create index to the recapture cohorts
 imatrix temp_data(1,nrecap,1,3); 
 //temp matrix to hold cohort index and nos. recaptures
 int ndx;
 for(int i = 1; i <= nrecap; i++)
 {
  int rel_size = rec_data(i,1);  

  for(int j = 1; j <= ncohort; j++)
   {
     if(rel_size == rel_data(j,2) )
     {
      ndx = rel_data(j,1); 
      break;
     }
   } //for j

   temp_data(i,1) = ndx;
   temp_data(i,2) = rec_data(i,1);
   temp_data(i,3) = rec_data(i,2);
 }  //for i 
 
 //stratify the recapture data 
 int ntime = max_days/bin_days + 1;  TRACE(ntime);
 imatrix obs_caps(1,ncohort,0,ntime);   //matrix for stratified data
 obs_caps.initialize();

 for(int rc = 1; rc <= nrecap; rc++)
 {
  int cc1 = temp_data(rc,1);  
  int da  = temp_data(rc,3);  
  int ti  = int( double(da)/double(bin_days) + 0.5 ); 
  if(ti > ntime)
  {
   ti = ntime; 
  }
  obs_caps(cc1,ti)++; 
   
 }
 TRACE(sum(obs_caps)); 

 //get min_size and max_size
 ivector temp(1,ncohort);
 temp = column(rel_data,2);
 int min_size = min(temp);
 int max_size = max(temp); 

 
 //output data
 ofstream out("yft.dat");
 out << "#min_size    max_size" << endl;
 out << min_size << "  " << max_size << endl;
 out <<"#ncohort" << endl;
 out << ncohort + 1<< endl;
 out <<"#release data" << endl;

 int sumrel = 0; 
 ivector sumrec(0,ntime);
 sumrec.initialize();

 for(int i = 1; i <= ncohort; i++)
 {
  sumrel += rel_data(i,3);
  out << rel_data(i,1) << "  "
      << rel_data(i,2) << "  " 
      << rel_data(i,3) << endl; 
      //<< rel_data(i,4) << endl;
 } 
 
 out << "#total release dummy cohort  " << endl;
 out << ncohort+1 <<  "  " << 76 << "  " << sumrel << endl;
  
 out << endl;
 out << "#recapture data" << endl;
 out << "#ntime"<< endl;
 out << ntime << endl;
 for(int i = 1; i <= ncohort; i++)
 {
  out << "#cohort " << i << endl;
  for(int j = 0; j <= ntime; j++)
  {
   sumrec(j) += obs_caps(i,j);
   out << "  " << obs_caps(i,j); 
  }
  out << endl;
 }
 out <<"#sum from all the cohorts " << endl; 
 out << sumrec  << endl;
   
 return 0; 
}
