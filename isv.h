/*
isv.h (31/1/03)
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2003, Darren Wilkinson
ISV: ANSI 'C' code library for stochastic volatility modelling

Header file for ISV code library
*/


/* includes */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>


/* main ISV object */
typedef struct
{
  int n;
  double phi;
  double mu;
  double sigma_eta;
  gsl_vector * y;
  gsl_vector * alpha_mean;
  gsl_vector * alpha_var;
  gsl_vector * alpha;
  gsl_vector * alpha_proposal;
} isv;


/* public common function prototypes */
isv * isv_alloc(double phi,double sigma_eta,double mu,gsl_vector *y);
void isv_filter(isv *self,int start,int finish);
void isv_simsmoother(isv *self,gsl_rng *r,int start,int finish);
double isv_loglik_norm(isv *self,gsl_vector *a,int start,int finish);
double isv_loglik(isv *self,gsl_vector *a,int start,int finish);
void isv_copy_alpha(isv *self,int start,int finish);
void isv_update_alpha(isv *self,gsl_rng *r,int start,int finish);
void isv_update_alpha_norm(isv *self,gsl_rng *r,int start,int finish);
double isv_sumdiff_alpha(isv *self,int start,int finish);
double isv_sumprod_alpha(isv *self,int start,int finish);
double isv_sumcrossprod_alpha(isv *self,int start,int finish);
double isv_ss_eta(isv *self,int start,int finish);
void isv_update_alpha_blocks(isv *self,gsl_rng *r,int blocksize);
void isv_update_mu(isv *self,gsl_rng *r,double a,double b);
void isv_update_phi(isv *self,gsl_rng *r,double a,double b);
void isv_update_phi_beta(isv *self,gsl_rng *r,double a,double b,double sd);
void isv_update_sigma_eta(isv *self,gsl_rng *r,double a,double b);
int isv_acceptlp(gsl_rng *r,double a);


/* public macros */
#define ISV_LETS_MEAN -1.27   /* more accurate would be better! */
#define ISV_LETS_VAR 4.934802201
#define ISV_PI 3.141592654
#define ISV_LTPOT 0.918938533 /* log(2 pi)/2 */


/* gsl macros */
#define GVG gsl_vector_get
#define GVS gsl_vector_set
#define GMG gsl_matrix_get
#define GMS gsl_matrix_set



/* eof */

