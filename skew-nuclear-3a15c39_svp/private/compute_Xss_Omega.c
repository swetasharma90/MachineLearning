#include "mex.h"
/* Computes b=P_Omega(UV'). 
   U is m\times k matrix.
   V is n\times k matrix.
   Omega is the vector of sampling indices.
   Usage: b=compute_X_Omega(U,V,Omega).
   Written by: Prateek Jain (pjain@cs.utexas.edu) and Raghu Meka (raghu@cs.utexas.edu)
   Last updated: 10/29/09
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

  mwSize m, n, k,Omega_size;
  mwIndex i,j,Omega_idx;
  double *U, *V, *Omega, *b;
  int l;
  m=mxGetM(prhs[0]);
  n=mxGetM(prhs[1]);
  k=mxGetN(prhs[0]);
  U=mxGetPr(prhs[0]);
  V=mxGetPr(prhs[1]);
  
  if (m != n) {
      mexErrMsgIdAndTxt("svt:skewSymmetricDims",
            "compute_Xss_Omega called with %i = size(U,1) != size(V,1) = %i",
             m, n);
  }

  Omega_size=mxGetM(prhs[2]);
  Omega=mxGetPr(prhs[2]);
  plhs[0]=mxCreateDoubleMatrix(Omega_size, 1, mxREAL);
  b=mxGetPr(plhs[0]);

  for(Omega_idx=0;Omega_idx<Omega_size;Omega_idx++){
    i=(mwSize)(((mwSize)Omega[Omega_idx]-1)%m);
    j=(mwSize)(((mwSize)Omega[Omega_idx]-1)/m);
    b[Omega_idx]=0;
    for(l=0;l<k;l++){
      b[Omega_idx]+=U[l*m+i]*V[l*n+j];
    }
    for(l=0;l<k;l++){
      b[Omega_idx]-=U[l*m+j]*V[l*n+i];
    }
    b[Omega_idx] *= 0.5;
  }
}
