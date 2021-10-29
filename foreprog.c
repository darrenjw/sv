/*
fsvprog.c
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2004, Darren Wilkinson
FSV: ANSI 'C' code library for stochastic volatility modelling

Example code for FSV models
*/

#include "forecast_fsv.h"

#define N 4

int main(int argc,char *argv[])
{
  gsl_matrix *y; 
  gsl_rng *r; int n,m,i; long it,iters;
  fsv *fsv; fcast *fcast; FILE *s;
  if (argc != 4) {
    fprintf(stderr,"Usage: %s <iters> <datapoints> <blocksize>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  iters=atoi(argv[1]); n=atoi(argv[2]); m=atoi(argv[3]);
  y=gsl_matrix_alloc(n,N);
  s=fopen("fsvdata.dat","r");
  if (s==NULL) {
    perror("failed to open data file");
    exit(EXIT_FAILURE);
  }
  gsl_matrix_fscanf(s,y);
  r=gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,655284);
  
  /*            phi sig  mu fphi fsig k y   */
  fsv=fsv_alloc(0.7,0.02,1.0,0.8,0.1,1,y);
  fcast=fcast_alloc(fsv,10); /* allow forecasts up to 10 steps */
  printf("Iter fphi fsig ");
  for (i=0;i<N;i++) {
    printf("phi[%d] sig[%d] mu[%d] b[%d] ",i,i,i,i);
  }
  for (i=0;i<N;i++) {
    printf("f[5,%d] ",i);
  }
  printf("\n");
  for (it=0;it<iters;it++) {
    fsv_update_f(fsv,r);
    fsv_update_alpha_f(fsv,r,m);
    /* fsv_update_phi_f(fsv,r,0.8,0.01); */
    fsv_update_phi_f_beta(fsv,r,1,1,0.05);
    fsv_update_sigma_f(fsv,r,1.0,0.01);
    fsv_update_alpha_w(fsv,r,m);
    fsv_update_phi_w(fsv,r,0.8,0.01);
    fsv_update_sigma_w(fsv,r,1.0,0.01);
    fsv_update_mu_w(fsv,r,0.0,100.0);
    fsv_update_b(fsv,r,1.0,100.0);
    fcast_sim(fcast,r); /* sample the forecasts (10 steps) */
    fcast_returns(fcast,5); /* look at 5 step returns */
    printf("%ld %f %f ",it,fsv->isvf[0]->phi,fsv->isvf[0]->sigma_eta);
    for (i=0;i<N;i++) {
      printf("%f %f %f %f ",fsv->isvw[i]->phi,fsv->isvw[i]->sigma_eta,
	     fsv->isvw[i]->mu,GMG(fsv->b,i,0));
    }
    for (i=0;i<N;i++) {
      printf("%f ",GVG(fcast->returns,i));
    }
    printf("\n");
  }
  return(EXIT_SUCCESS);
}

/* eof */

