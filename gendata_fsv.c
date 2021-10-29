/*
gendata_fsv.c
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2004, Darren Wilkinson
ISV: ANSI 'C' code library for stochastic volatility modelling

Generate a simulated data set based on a FSV model
*/

#include "fsv.h"


#define N 4

/* function prototypes */

void sim_isv(gsl_rng *r,gsl_vector *y,double phi,double sigma,double mu);


/* main function */

int main(int argc,char *argv[])
{
  gsl_vector *f,*b;
  gsl_vector_view v;
  gsl_matrix *y;
  gsl_rng *r; 
  int i,t,n;
  if (argc != 2) {
    fprintf(stderr,"Usage: %s <datapoints>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  n=atoi(argv[1]);

  r=gsl_rng_alloc(gsl_rng_mt19937);
  f=gsl_vector_alloc(n);
  b=gsl_vector_alloc(N);
  y=gsl_matrix_alloc(n,N);

  sim_isv(r,f,0.8,0.1,0.0);
  for (i=0;i<N;i++) {
    v=gsl_matrix_column(y,i);
    sim_isv(r,&(v.vector),0.7,0.02,i);
    GVS(b,i,1.0/(i+1.0));
  }
  gsl_blas_dger(1.0,f,b,y);
  for (t=0;t<n;t++) {
    for (i=0;i<N;i++) {
      printf("%f ",GMG(y,t,i));
    }
    printf("\n");
  }
  return(EXIT_SUCCESS);
}


/* helper functions */

void sim_isv(gsl_rng *r,gsl_vector *y,double phi,double sigma,double mu)
{
  int i;
  double v,alpha;
  v=sigma*sigma/(1.0-phi*phi);
  alpha=mu+gsl_ran_gaussian(r,sqrt(v));
  for (i=0;i<y->size;i++) {
    alpha=mu+phi*(alpha-mu)+gsl_ran_gaussian(r,sigma);
    GVS(y,i,gsl_ran_gaussian(r,exp(0.5*alpha)));
  }
}


/* eof */

