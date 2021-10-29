/*
isvprog.c
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2003, Darren Wilkinson
ISV: ANSI 'C' code library for stochastic volatility modelling

Example code for ISV models
*/

#include "isv.h"

int main(int argc,char *argv[])
{
  gsl_vector *y; 
  gsl_rng *r; int n,m,i; long it,iters;
  isv *isv; FILE *s;
  if (argc != 4) {
    fprintf(stderr,"Usage: %s <iters> <datapoints> <blocksize>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  iters=atoi(argv[1]); n=atoi(argv[2]); m=atoi(argv[3]);
  y=gsl_vector_alloc(n);
  s=fopen("isvdata.dat","r");
  if (s==NULL) {
    perror("failed to open data file");
    exit(EXIT_FAILURE);
  }
  gsl_vector_fscanf(s,y);
  r=gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,655284);
  
  /*            phi sig  mu y      */
  isv=isv_alloc(0.8,0.1,1.0,y);
  /* isv=isv_alloc(0.5,1.5,2.0,y); */
  isv_update_alpha_norm(isv,r,0,n-1);
  printf("Iter mu phi sigma_eta");
  /* for (i=0;i<n;i++) { */
  for (i=0;i<10;i++) {
    printf(" alpha[%d]",i);
  }
  printf("\n");
  for (it=0;it<iters;it++) {
    isv_update_alpha_blocks(isv,r,m);
    isv_update_mu(isv,r,0,10000);
    isv_update_phi(isv,r,0.8,0.01);
    /* isv_update_sigma_eta(isv,r,1.0,0.001); */
    isv_update_sigma_eta(isv,r,0.001,0.001);
    printf("%ld %f %f %f",it,
	   isv->mu,isv->phi,isv->sigma_eta);
    /* for (i=0;i<n;i++) { */
    for (i=0;i<10;i++) {
      printf(" %f",GVG(isv->alpha,i));
    }
    printf("\n");
  }

  return(EXIT_SUCCESS);
}

/* eof */

