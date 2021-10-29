/*
gendata_isv.c
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2003, Darren Wilkinson
ISV: ANSI 'C' code library for stochastic volatility modelling

Generate a simulated data set based on an ISV model
*/

#include "isv.h"

int main(int argc,char *argv[])
{
  double m,v,true_phi,true_sigma_eta,true_mu;
  gsl_vector *true_alpha,*y; 
  gsl_rng *r; int i,n;
  if (argc != 2) {
    fprintf(stderr,"Usage: %s <datapoints>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  n=atoi(argv[1]);
  true_phi=0.8; true_sigma_eta=0.1; true_mu=1.0;
  true_alpha=gsl_vector_alloc(n);
  y=gsl_vector_alloc(n);
  r=gsl_rng_alloc(gsl_rng_mt19937);
  m=true_mu;
  v=true_sigma_eta*true_sigma_eta/(1-true_phi*true_phi);
  GVS(true_alpha,0,m+gsl_ran_gaussian(r,sqrt(v)));
  GVS(y,0,gsl_ran_gaussian(r,exp(gsl_vector_get(true_alpha,0)/2)));
  for (i=1;i<n;i++) {
    m=true_mu+true_phi*(gsl_vector_get(true_alpha,i-1)-true_mu);
    GVS(true_alpha,i,m+gsl_ran_gaussian(r,true_sigma_eta));
    GVS(y,i,gsl_ran_gaussian(r,exp(GVG(true_alpha,i)/2)));
  }
  fprintf(stderr,"true_sigma_eta:%f, true_mu:%f, true_phi:%f,\ntrue_alpha[0]:%f, true_alpha[100]:%f\n",true_sigma_eta,true_mu,true_phi,GVG(true_alpha,0),GVG(true_alpha,100));
  
  gsl_vector_fprintf(stdout,y,"%f");

  return(EXIT_SUCCESS);
}

/* eof */

