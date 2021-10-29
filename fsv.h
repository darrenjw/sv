/*
fsv.h (18/5/04)
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2004, Darren Wilkinson
FSV: ANSI 'C' code library for stochastic volatility modelling

Header file for FSV code library
*/


#include "isv.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

/* main FSV object */
typedef struct
{
  int k, n, nn;
  gsl_matrix *b, *f, *y;
  isv **isvf, **isvw;
} fsv;


/* public common function prototypes */
fsv * fsv_alloc(double phi, double sig, double mu, double fphi,
		double fsig, int k, gsl_matrix *y);
void fsv_update_alpha_f(fsv *self,gsl_rng *r,int blocksize);
void fsv_update_alpha_w(fsv *self,gsl_rng *r,int blocksize);
void fsv_update_sigma_f(fsv *self,gsl_rng *r,double a,double b);
void fsv_update_phi_f(fsv *self,gsl_rng *r,double a,double b);
void fsv_update_phi_f_beta(fsv *self,gsl_rng *r,double a,double b,double sd);
void fsv_update_sigma_w(fsv *self,gsl_rng *r,double a,double b);
void fsv_update_phi_w(fsv *self,gsl_rng *r,double a,double b);
void fsv_update_phi_w_beta(fsv *self,gsl_rng *r,double a,double b,double sd);
void fsv_update_mu_w(fsv *self,gsl_rng *r,double a,double b);
void fsv_canorm_inplace(gsl_rng *r,gsl_vector *h,gsl_matrix *k);
double fsv_canorm_scalar(gsl_rng *r,double h,double k);
void fsv_update_f(fsv *self,gsl_rng *r);
void fsv_update_b(fsv *self,gsl_rng *r,double a,double b);


/* eof */

