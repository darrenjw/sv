/*
forecast_isv.h
Author: Darren Wilkinson <d.j.wilkinson@ncl.ac.uk>
<URL:http://www.staff.ncl.ac.uk/d.j.wilkinson/>
Copyright (C) 2004, Darren Wilkinson
FSV: ANSI 'C' code library for stochastic volatility modelling

Headers for forecasting code

*/

#include "fsv.h"

/* forecast object */
typedef struct
{
  fsv *fsv;
  gsl_matrix *foremat,*fmat;
  gsl_vector *returns;
  int k,nn,T;
} fcast;

/* public function prototypes */

fcast * fcast_alloc(fsv *fsv,int T);
void fcast_sim(fcast *self,gsl_rng *r);
void fcast_returns(fcast *self,int t);


/* eof */

