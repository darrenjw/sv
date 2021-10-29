/*
test_isv.c
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2003, Darren Wilkinson
ISV: ANSI 'C' code library for stochastic volatility modelling

Test file for ISV code library
*/

#include "isv.h"

int dtest(char * name,double should,double is,double tol);

int main(void)
{
  isv *isv; gsl_vector *data; gsl_rng *r;
  double phi,sigma_eta,mu,sum,ss,val;
  long i,j,n; int err;
  err=0; n=1e6;
  phi=0.8; sigma_eta=0.1; mu=1;
  data=gsl_vector_alloc(5);
  GVS(data,0,5);
  GVS(data,1,-4);
  GVS(data,2,3);
  GVS(data,3,1);
  GVS(data,4,-1);
  isv=isv_alloc(phi,sigma_eta,mu,data);
  isv_filter(isv,0,4);
  err+=dtest("filter mean",1.027748,GVG(isv->alpha_mean,4),1e-5);
  err+=dtest("filter var",0.0273973,GVG(isv->alpha_var,4),1e-5);
  
  sum=0; ss=0;
  r=gsl_rng_alloc(gsl_rng_mt19937);
  for (i=0;i<n;i++) {
    isv_simsmoother(isv,r,0,4);
    val=GVG(isv->alpha_proposal,0);
    sum+=val;
    ss+=(val-1.0429)*(val-1.0429);
  }
  err+=dtest("simsmoother mean",1.0429,sum/n,1e-3);
  err+=dtest("simsmoother var",0.0274,ss/n,1e-4);

  sum=0; ss=0;
  for (i=0;i<n;i++) {
    isv_filter(isv,0,2);
    isv_simsmoother(isv,r,0,2);
    isv_copy_alpha(isv,0,2);
    isv_filter(isv,0,2);
    isv_simsmoother(isv,r,3,4);
    isv_copy_alpha(isv,3,4);
    val=GVG(isv->alpha,0);
    sum+=val;
    ss+=(val-1.0429)*(val-1.0429);
  }
  err+=dtest("block update mean",1.0429,sum/n,1e-3);
  err+=dtest("block update var",0.0274,ss/n,1e-4);

  sum=0; ss=0;
  for (i=0;i<n;i++) {
    for (j=0;j<5;j++) {
      isv_filter(isv,j,j);
      isv_simsmoother(isv,r,j,j);
      isv_copy_alpha(isv,j,j);
    }
    val=GVG(isv->alpha,0);
    sum+=val;
    ss+=(val-1.0429)*(val-1.0429);
  }
  err+=dtest("sequential update mean",1.0429,sum/n,1e-3);
  err+=dtest("sequential update var",0.0274,ss/n,1e-4);

  sum=0; ss=0;
  for (i=0;i<n;i++) {
    for (j=0;j<5;j+=2) {
      isv_filter(isv,j,j);
      isv_simsmoother(isv,r,j,j);
      isv_copy_alpha(isv,j,j);
    }
    for (j=1;j<5;j+=2) {
      isv_filter(isv,j,j);
      isv_simsmoother(isv,r,j,j);
      isv_copy_alpha(isv,j,j);
    }
    val=GVG(isv->alpha,0);
    sum+=val;
    ss+=(val-1.0429)*(val-1.0429);
  }
  err+=dtest("odd-even update mean",1.0429,sum/n,1e-3);
  err+=dtest("odd-even update var",0.0274,ss/n,1e-4);


   
  if (err==0) {
    printf("All tests passed.\n");
    return(EXIT_SUCCESS);
  } else {
    fprintf(stderr,"ERROR: %d tests failed\n",err);
    printf("ERROR: %d tests failed\n",err);
    return(EXIT_FAILURE);
  }
}

int dtest(char * name,double should,double is,double tol)
{
  int err;
  if ( fabs(should-is) > tol ) {
    fprintf(stderr,"FAILED: Test: %s, Should: %f, Is: %f, Tol: %f\n",
	    name,should,is,tol);
    printf("FAILED: Test: %s, Should: %f, Is: %f, Tol: %f\n",
	    name,should,is,tol);
    err=1;
  } else {
    printf("Test :%s: passed.\n",name);
    err=0;
  }
  return(err);
}


/* eof */

