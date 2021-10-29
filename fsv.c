/*
fsv.c
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2004, Darren Wilkinson
FSV: ANSI 'C' code library for stochastic volatility modelling

Main code file for FSV library

*/

#include "fsv.h"


/* common function definitions */

fsv * fsv_alloc(double phi, double sig, double mu, double fphi,
		double fsig, int k, gsl_matrix *y)
{
  int i;
  fsv *self;
  self=malloc(sizeof(fsv));
  /* TODO: should check pointers! */
  /* dimensions */
  self->n=y->size1;
  self->nn=y->size2;
  self->k=k;
  /* set up factor vol series */
  self->isvf=malloc(k*sizeof(isv *));
  for (i=0;i<k;i++) {
    self->isvf[i]=isv_alloc(fphi,fsig,0.0,gsl_vector_calloc(self->n));
  }
  /* set up the ideosynchratic volatility series */
  self->isvw=malloc(self->nn*sizeof(isv *));
  for (i=0;i<self->nn;i++) {
    self->isvw[i]=isv_alloc(phi,sig,mu,gsl_vector_calloc(self->n));
  }
  /* other inits */
  self->y=y;
  self->f=gsl_matrix_calloc(self->n,k);
  self->b=gsl_matrix_calloc(self->nn,k);
  gsl_matrix_set_identity(self->b);
  return(self);
}

void fsv_update_alpha_f(fsv *self,gsl_rng *r,int blocksize)
{
  int i,j;
  for (i=0;i<self->k;i++) {
    for (j=0;j<self->n;j++) {
      GVS(self->isvf[i]->y,j,GMG(self->f,j,i));
    }
    isv_update_alpha_blocks(self->isvf[i],r,blocksize);
  }
}

void fsv_update_alpha_w(fsv *self,gsl_rng *r,int blocksize)
{
  gsl_matrix *work;
  int t,j;
  /* TODO: should make work an FSV attr */
  work=gsl_matrix_alloc(self->n,self->nn);
  /* work=f*b' */
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,self->f,self->b,0.0,work);
  for (t=0;t<self->n;t++) {
    for (j=0;j<self->nn;j++) {
      GVS(self->isvw[j]->y,t,GMG(self->y,t,j)-GMG(work,t,j));
    }
  }
  gsl_matrix_free(work);
  for (j=0;j<self->nn;j++) {
    isv_update_alpha_blocks(self->isvw[j],r,blocksize);
  }
}

void fsv_update_sigma_f(fsv *self,gsl_rng *r,double a,double b)
{
  int i;
  for (i=0;i<self->k;i++) {
    isv_update_sigma_eta(self->isvf[i],r,a,b);
  }
}

void fsv_update_phi_f(fsv *self,gsl_rng *r,double a,double b)
{
  int i;
  for (i=0;i<self->k;i++) {
    isv_update_phi(self->isvf[i],r,a,b);
  }
}

void fsv_update_phi_f_beta(fsv *self,gsl_rng *r,double a,double b,double sd)
{
  int i;
  for (i=0;i<self->k;i++) {
    isv_update_phi_beta(self->isvf[i],r,a,b,sd);
  }
}

void fsv_update_sigma_w(fsv *self,gsl_rng *r,double a,double b)
{
  int j;
  for (j=0;j<self->nn;j++) {
    isv_update_sigma_eta(self->isvw[j],r,a,b);
  }
}

void fsv_update_mu_w(fsv *self,gsl_rng *r,double a,double b)
{
  int j;
  for (j=0;j<self->nn;j++) {
    isv_update_mu(self->isvw[j],r,a,b);
  }
}

void fsv_update_phi_w(fsv *self,gsl_rng *r,double a,double b)
{
  int j;
  for (j=0;j<self->nn;j++) {
    isv_update_phi(self->isvw[j],r,a,b);
  }
}

void fsv_update_phi_w_beta(fsv *self,gsl_rng *r,double a,double b,double sd)
{
  int j;
  for (j=0;j<self->nn;j++) {
    isv_update_phi_beta(self->isvw[j],r,a,b,sd);
  }
}

void fsv_canorm_inplace(gsl_rng *r,gsl_vector *h,gsl_matrix *k)
{
  int i;
  gsl_linalg_cholesky_decomp(k);
  gsl_blas_dtrsv(CblasLower,CblasNoTrans,CblasNonUnit,k,h);
  for (i=0;i<h->size;i++) {
    GVS(h,i,GVG(h,i)+gsl_ran_ugaussian(r));
  }
  gsl_blas_dtrsv(CblasLower,CblasTrans,CblasNonUnit,k,h);
}

double fsv_canorm_scalar(gsl_rng *r,double h,double k)
{
  return( h/k+gsl_ran_gaussian(r,1.0/sqrt(k)) );
}

void fsv_update_f(fsv *self,gsl_rng *r)
{
  gsl_matrix *tb, *fmat;
  gsl_vector_view view,view2;
  int i,t;
  /* TODO: should make tb and fmat FSV attrs */
  tb=gsl_matrix_alloc(self->nn,self->k);
  fmat=gsl_matrix_alloc(self->k,self->k);
  for (t=0;t<self->n;t++) {
    /* canonical location */
    gsl_matrix_memcpy(tb,self->b);
    for (i=0;i<self->nn;i++) {
      view=gsl_matrix_row(tb,i);
      gsl_vector_scale(&(view.vector),exp(-GVG(self->isvw[i]->alpha,t)));
    }
    view=gsl_matrix_row(self->f,t);
    view2=gsl_matrix_row(self->y,t);
    gsl_blas_dgemv(CblasTrans,1.0,tb,&(view2.vector),0.0,&(view.vector));
    /* precision */
    gsl_matrix_set_zero(fmat);
    for (i=0;i<self->k;i++) {
      GMS(fmat,i,i,exp(-GVG(self->isvf[i]->alpha,t)));
    }
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,self->b,tb,1.0,fmat);
    fsv_canorm_inplace(r,&(view.vector),fmat);
  }
  gsl_matrix_free(tb);
  gsl_matrix_free(fmat);
}
  
void fsv_update_b(fsv *self,gsl_rng *r,double a,double b)
{
  int l,m,mmax,t,i;
  double term,h,k;
  for (l=1;l<self->nn;l++) {
    mmax=l;
    if (mmax > self->k) mmax=self->k;
    for (m=0;m<mmax;m++) {
      h=a/b;
      k=1.0/b;
      for (t=0;t<self->n;t++) {
	term=0;
	for (i=0;i<k;i++) {
	  if (i==m) break;
	  term+=GMG(self->f,t,i)*GMG(self->b,l,i);
	}
	h+=GMG(self->f,t,m)*exp(-GVG(self->isvw[l]->alpha,t))
	  *(GMG(self->y,t,l)-term);
	k+=GMG(self->f,t,m)*GMG(self->f,t,m)*exp(-GVG(self->isvw[l]->alpha,t));
      }
      GMS(self->b,l,m,fsv_canorm_scalar(r,h,k));
    }
  }
}
		    


/* eof */

