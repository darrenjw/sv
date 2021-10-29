/*
isv.c
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2004, Darren Wilkinson
ISV: ANSI 'C' code library for stochastic volatility modelling

Main code file for ISV library

*/

#include "isv.h"


/* common function definitions */

isv * isv_alloc(double phi,double sigma_eta,double mu,gsl_vector *y)
{
  isv *self;
  self=malloc(sizeof(isv));
  /* TODO: should check pointers! */
  self->phi = phi;
  self->sigma_eta = sigma_eta;
  self->mu = mu;
  self->y = y;
  self->n = y->size;
  self->alpha_mean = gsl_vector_calloc(self->n);
  self->alpha_var = gsl_vector_calloc(self->n);
  self->alpha = gsl_vector_calloc(self->n);
  self->alpha_proposal = gsl_vector_alloc(self->n);
  return(self);
}

void isv_filter(isv *self,int start,int finish)
{
  int i;
  double m,v,yp;
#ifdef ISV_CHK
  if (!(
	(start >= 0) && (finish >= start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_filter pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  if (start==0) {
    m=(self->mu);
    v=(self->sigma_eta)*(self->sigma_eta)
      / (1.0 - (self->phi)*(self->phi));
  } else {
    m=(self->mu)
      + (self->phi)*(GVG(self->alpha,start-1)-(self->mu));
    v=(self->sigma_eta)*(self->sigma_eta);
  }
  for (i=start;i<=finish;i++) {
    yp=log(GVG(self->y,i)*GVG(self->y,i));
    m=(ISV_LETS_VAR*m + v*(yp-ISV_LETS_MEAN)) /
      (v+ISV_LETS_VAR);
    v=v*ISV_LETS_VAR/(v+ISV_LETS_VAR);
    GVS(self->alpha_mean,i,m);
    GVS(self->alpha_var,i,v);
    m=(self->mu)+(self->phi)*(m-(self->mu));
    v=(self->phi)*(self->phi)*v + (self->sigma_eta)*(self->sigma_eta);
  }
}

void isv_simsmoother(isv *self,gsl_rng *r,int start,int finish)
{
  int i;
  double m,v;
#ifdef ISV_CHK
  if (!(
	(start >= 0) && (finish >= start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_simsmoother pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  if (finish == (self->n)-1) {
    GVS(self->alpha_proposal,finish,
		   GVG(self->alpha_mean,finish) + 
		   gsl_ran_gaussian(r,sqrt(GVG(self->alpha_var,finish))));
  } else {
    m=GVG(self->alpha_mean,finish) +
      (self->phi)*GVG(self->alpha_var,finish)*
      (GVG(self->alpha,finish+1)-(self->mu)-(self->phi)*
       (GVG(self->alpha_mean,finish)-(self->mu)))/
      ((self->phi)*(self->phi)*GVG(self->alpha_var,finish)+
       (self->sigma_eta)*(self->sigma_eta));
    v=GVG(self->alpha_var,finish)*(self->sigma_eta)*(self->sigma_eta)/
      ((self->phi)*(self->phi)*GVG(self->alpha_var,finish)+
       (self->sigma_eta)*(self->sigma_eta));
    GVS(self->alpha_proposal,finish,m+gsl_ran_gaussian(r,sqrt(v)));    
  }
  for (i=finish-1;i>=start;i--) {
    m=GVG(self->alpha_mean,i) +
      (self->phi)*GVG(self->alpha_var,i)*
      (GVG(self->alpha_proposal,i+1)-(self->mu)-(self->phi)*
       (GVG(self->alpha_mean,i)-(self->mu)))/
      ((self->phi)*(self->phi)*GVG(self->alpha_var,i)+
       (self->sigma_eta)*(self->sigma_eta));
    v=GVG(self->alpha_var,i)*(self->sigma_eta)*(self->sigma_eta)/
      ((self->phi)*(self->phi)*GVG(self->alpha_var,i)+
       (self->sigma_eta)*(self->sigma_eta));
    GVS(self->alpha_proposal,i,m+gsl_ran_gaussian(r,sqrt(v)));
  }
}

double isv_loglik_norm(isv *self,gsl_vector *a,int start,int finish)
{
  double sum,term,yp;
  int i;
#ifdef ISV_CHK
  if (!(
	(start >= 0) && (finish >= start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_loglik_norm pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  sum=0;
  for (i=start;i<=finish;i++) {
    yp=log(GVG(self->y,i)*GVG(self->y,i));
    term=yp-GVG(a,i)-ISV_LETS_MEAN;
    sum+=term*term;
  }
  return( -1.5*(1+finish-start)*log(ISV_PI) - sum/(ISV_PI*ISV_PI) );
}

double isv_loglik(isv *self,gsl_vector *a,int start,int finish)
{
  double sum,yi;
  int i;
#ifdef ISV_CHK
  if (!(
	(start >= 0) && (finish >= start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_loglik pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  sum=0;
  for (i=start;i<=finish;i++) {
    yi=GVG(self->y,i);
    sum+=GVG(a,i) + yi*yi/exp(GVG(a,i));
  }
  return( ISV_LTPOT*(1+finish-start) - 0.5*sum );
}

void isv_copy_alpha(isv *self,int start,int finish)
{
  int i;
#ifdef ISV_CHK
  if (!(
	(start >= 0) && (finish >= start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_copy_alpha pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  for (i=start;i<=finish;i++) {
    GVS(self->alpha,i, GVG(self->alpha_proposal,i) );
  }
}

void isv_update_alpha(isv *self,gsl_rng *r,int start,int finish)
{
  double a;
#ifdef ISV_CHK
  if (!(
	(start >= 0) && (finish >= start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_update_alpha pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  isv_filter(self,start,finish);
  isv_simsmoother(self,r,start,finish);
  a = isv_loglik(self,self->alpha_proposal,start,finish)
    -isv_loglik_norm(self,self->alpha_proposal,start,finish)
    -(isv_loglik(self,self->alpha,start,finish)
      -isv_loglik_norm(self,self->alpha,start,finish));
  if ( isv_acceptlp(r,a) )
    isv_copy_alpha(self,start,finish);
}

void isv_update_alpha_norm(isv *self,gsl_rng *r,int start,int finish)
{
#ifdef ISV_CHK
  if (!(
	(start >= 0) && (finish > start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_update_alpha pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  isv_filter(self,start,finish);
  isv_simsmoother(self,r,start,finish);
  isv_copy_alpha(self,start,finish);
}

double isv_sumdiff_alpha(isv *self,int start,int finish)
{
  double sum; int i;
#ifdef ISV_CHK
  if (!(
	(start > 0) && (finish >= start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_sumdiff_alpha pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  sum=0;
  for (i=start;i<=finish;i++) {
    sum+=(GVG(self->alpha,i)
	  -(self->phi)*GVG(self->alpha,i-1));
  }
  return(sum);
}

double isv_sumprod_alpha(isv *self,int start,int finish)
{
  double sum; int i;
#ifdef ISV_CHK
  if (!(
	(start > 0) && (finish > start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_sumprod_alpha pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  sum=0;
  for (i=start;i<=finish;i++) {
    sum+=(GVG(self->alpha,i)-(self->mu))*
      (GVG(self->alpha,i)-(self->mu));
  }
  return(sum);
}

double isv_sumcrossprod_alpha(isv *self,int start,int finish)
{
  double sum; int i;
#ifdef ISV_CHK
  if (!(
	(start > 0) && (finish >= start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_sumcrossprod_alpha pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  sum=0;
  for (i=start;i<=finish;i++) {
    sum+=(GVG(self->alpha,i)-(self->mu))*
      (GVG(self->alpha,i-1)-(self->mu));
  }
  return(sum);
}

double isv_ss_eta(isv *self,int start,int finish)
{
  double sum,term; int i;
#ifdef ISV_CHK
  if (!(
	(start > 0) && (finish > start) && (finish < self->n)
	)) {
    fprintf(stderr,"isv_ss_eta pre-conditions failed\n");
    exit(EXIT_FAILURE);
  }
#endif
  sum=0;
  for (i=start;i<=finish;i++) {
    term=(GVG(self->alpha,i)-(self->mu)) 
      - (self->phi)*(GVG(self->alpha,i-1)-(self->mu));
    sum+=term*term;
  }
  /* fprintf(stderr,"start=%d, finish=%d, mu=%f, phi=%f, sum=%f\n",start,finish,self->mu,self->phi,sum); */
  return(sum);
}


void isv_update_alpha_blocks_alt(isv *self,gsl_rng *r,int blocksize)
{
  int start;
  start=0;
  while (start<(self->n)) {
    isv_update_alpha(self,r,start,start+blocksize-1);
    start+=2*blocksize;
  }
  start=blocksize;
  while (start<(self->n)) {
    isv_update_alpha(self,r,start,start+blocksize-1);
    start+=2*blocksize;
  }
}

void isv_update_alpha_blocks_seq(isv *self,gsl_rng *r,int blocksize)
{
  int start;
  start=0;
  while (start<(self->n)) {
    isv_update_alpha(self,r,start,start+blocksize-1);
    start+=blocksize;
  }
}

void isv_update_alpha_blocks(isv *self,gsl_rng *r,int blocksize)
{
  int start;
  start=0;
  while (start<(self->n)) {
    if ( start + blocksize > (self->n) ) {
      isv_update_alpha(self,r,start,(self->n)-1);
    } else {
      isv_update_alpha(self,r,start,start+blocksize-1);
    }
    start+=blocksize;
  }
}

void isv_update_mu(isv *self,gsl_rng *r,double a,double b)
{
  /* Assuming a N(a,b) prior on mu */
  double m,v,p,sumdiff;
  sumdiff=isv_sumdiff_alpha(self,1,(self->n)-1);
  p=(1.0/b) + ((self->n)-1)*(1-(self->phi))*(1-(self->phi))/
    ((self->sigma_eta)*(self->sigma_eta)) ;
  v=1.0/p;
  m=v*( (a/b) + (1-(self->phi))*sumdiff/
	((self->sigma_eta)*(self->sigma_eta)) );
  (self->mu)=m+gsl_ran_gaussian(r,sqrt(v));
}

void isv_update_phi(isv *self,gsl_rng *r,double a,double b)
{
  /* Assuming a N(a,b) prior on phi */
  double m,v,p,sump,sumc;
  sump=isv_sumprod_alpha(self,1,(self->n)-1);
  sumc=isv_sumcrossprod_alpha(self,1,(self->n)-1);
  p=(1.0/b) + sump/((self->sigma_eta)*(self->sigma_eta));
  v=1.0/p;
  m=v*( (a/b) + sumc/((self->sigma_eta)*(self->sigma_eta)) );
  (self->phi)=m+gsl_ran_gaussian(r,sqrt(v));
}

void isv_update_phi_beta(isv *self,gsl_rng *r,double a,double b,double sd)
{
  /* Assuming a Beta(a,b) prior on phi and a random walk 
     Metropolis-Hastings update */
  double m,v,aprob,sump,sumc,phistar;
  phistar=(self->phi)+gsl_ran_gaussian(r,sd);
  if ( (phistar<0.0) || (phistar>1.0) )
    return;
  sump=isv_sumprod_alpha(self,1,(self->n)-1);
  sumc=isv_sumcrossprod_alpha(self,1,(self->n)-1);
  v=sump/((self->sigma_eta)*(self->sigma_eta));
  m=sumc/sump;
  aprob=(-0.5*v)*(((phistar-m)*(phistar-m))
		  -((self->phi - m)*(self->phi - m)))
    + (a-1)*((log(phistar)) - log(self->phi)) 
    + (b-1)*((log(1-phistar)) - log(1-self->phi));
  if( isv_acceptlp(r,aprob) )
    (self->phi) = phistar;
}

void isv_update_sigma_eta(isv *self,gsl_rng *r,double a,double b)
{
  /* Under a Gamma(a,b) prior on the corresponding precision */
  double p,ss;
  ss=isv_ss_eta(self,1,(self->n)-1);
  p=gsl_ran_gamma(r,a+0.5*((self->n)-1),1.0/(b+0.5*ss));
  /* fprintf(stderr,"n=%d, a=%f, b=%f, ss=%f, p=%f\n",self->n,a,b,ss,p); */
  (self->sigma_eta)=1.0/sqrt(p);
}

int isv_acceptlp(gsl_rng *r,double a)
{
  return( log(gsl_rng_uniform(r)) < a );
}

/* eof */

