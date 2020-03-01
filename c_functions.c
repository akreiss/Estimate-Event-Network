#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>


SEXP intensities(SEXP theta,SEXP Y)
{
  // Initialize Data
  int dim,k,N,r;
  SEXP intensities;
  double thetaY;
  
  dim=LENGTH(theta);
  N=LENGTH(Y);
  
  intensities=PROTECT(NEW_NUMERIC(N));
  
  // Go through all covariates and compute the intensity
  for(k=0;k<=N-1;k++)
  {
    thetaY=0;
    for(r=0;r<=dim-1;r++)
      thetaY=thetaY+REAL(theta)[r]*REAL(VECTOR_ELT(Y,k))[r];
    
    REAL(intensities)[k]=exp(thetaY);
  }
  
  UNPROTECT(1);
  
  return(intensities);
}