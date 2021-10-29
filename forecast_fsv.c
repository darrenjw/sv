/*
forecast_isv.c
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2004, Darren Wilkinson
FSV: ANSI 'C' code library for stochastic volatility modelling

Forecasting code

*/

#include "forecast_fsv.h"

/* private function prototypes */
void fcast_sim(fcast *self,gsl_rng *r);
void sim_isv(gsl_rng *r,gsl_vector *y,
	     double phi,double sigma,double mu,double alpha);


/* function declarations */

fcast * fcast_alloc(fsv *fsv,int T)
{
  fcast *self;
  self=malloc(sizeof(fcast));
  self->foremat=gsl_matrix_alloc(T,fsv->nn);
  self->fmat=gsl_matrix_alloc(T,fsv->k);
  self->returns=gsl_vector_alloc(fsv->nn);
  /* TODO: should check for null pointers! */
  self->fsv=fsv;
  self->T=T;
  self->nn=fsv->nn;
  self->k=fsv->k;
  return(self);
}

void fcast_returns(fcast *self,int t)
{
  /* assumes that one-step returns have already been simulated */
  int i,j;
  if (t>self->T) {
    fprintf(stderr,"Can't forecast returns more than %d time points ahead\n",self->T);
    exit(EXIT_FAILURE);
  }
  gsl_vector_set_zero(self->returns);
  for (j=0;j<t;j++) {
    for (i=0;i<self->nn;i++) {
      GVS(self->returns,i,GVG(self->returns,i)+GMG(self->foremat,j,i));
    }
  }
}

void fcast_sim(fcast *self,gsl_rng *r)
{
  int i,j;
  gsl_vector_view v;
  isv *isv;
  for (j=0;j<self->k;j++) {
    isv=self->fsv->isvf[j];
    v=gsl_matrix_column(self->fmat,j);
    sim_isv(r,&(v.vector),isv->phi,isv->sigma_eta,isv->mu,
	    GVG(isv->alpha,(self->nn)-1));
  }
  for (i=0;j<self->nn;j++) {
    isv=self->fsv->isvw[j];
    v=gsl_matrix_column(self->foremat,j);
    sim_isv(r,&(v.vector),isv->phi,isv->sigma_eta,isv->mu,
	    GVG(isv->alpha,(self->nn)-1));
  }
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,self->fmat,
		 self->fsv->b,1.0,self->foremat);
}


void sim_isv(gsl_rng *r,gsl_vector *y,
	     double phi,double sigma,double mu,double alpha)
{
  int i;
  for (i=0;i<y->size;i++) {
    alpha=mu+phi*(alpha-mu)+gsl_ran_gaussian(r,sigma);
    GVS(y,i,gsl_ran_gaussian(r,exp(0.5*alpha)));
  }
}




/* eof  */

