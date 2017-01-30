/**
 * Compute a skew-symmetric sparse matrix from a set of ratings.
 * This function produces either a dense matrix.
 * @author David F. Gleich
 *
 * History
 * -------
 * :2009-11-13: Initial coding
 * :2009-12-05: Added log-odds function
 * :2009-12-06: Fixed transpose from log-odds function.
 */

#include <mex.h>
#include "mexhelp.h"
#include <math.h>

typedef struct {
    mwIndex* cp;
    mwIndex* ri;
    double* a;
    mwSize m;
    mwSize n;
} cs;

void cs_mex_get_sparse(const mxArray* Amatlab, cs *A) {
    A->cp = mxGetJc(Amatlab);
    A->ri = mxGetIr(Amatlab);
    A->a = mxGetPr(Amatlab);
    A->m = mxGetM(Amatlab);
    A->n = mxGetN(Amatlab);
}

void skew_arithmean(cs *At, double* Y, double* N) {
    mwSize nusers = At->n;
    mwSize nitems = At->m;
    mwIndex i, j, k, mi, mj;
    double ri, rj;
    /* allocate memory for logs, they are really slow unless precomputed */
    mwIndex *ms = calloc(nitems, sizeof(mwIndex));
    double *rs = calloc(nitems, sizeof(double));
    for (i=0; i<nusers; i++) {
        /* examine column i of At, which has the ratings for user i */
        mwIndex nrats = At->cp[i+1]-At->cp[i];
        for (j=0; j<nrats; j++) {
            rs[j] = At->a[At->cp[i]+j];
            ms[j] = At->ri[At->cp[i]+j];
        }
        for (j=0; j<nrats; j++) {
            mi=ms[j];
            ri=rs[j];
            for (k=0; k<nrats; k++) {
                if (k==j) { continue; } /* skip the same rating */
                mj=ms[k];
                rj=rs[k];
                Y[nitems*mj+mi] += ri-rj;
                N[nitems*mj+mi]++;
            }
        }
    }
    free(ms);
    free(rs);
}

void skew_geomean(cs *At, double* Y, double* N, double shift) {
    mwSize nusers = At->n;
    mwSize nitems = At->m;
    mwIndex i, j, k, mi, mj;
    double ri, rj;
    /* allocate memory for logs, they are really slow unless precomputed */
    mwIndex *ms = calloc(nitems, sizeof(mwIndex));
    double *rs = calloc(nitems, sizeof(double));
    for (i=0; i<nusers; i++) {
        /* examine column i of At, which has the ratings for user i */
        mwIndex nrats = At->cp[i+1]-At->cp[i];
        for (j=0; j<nrats; j++) {
            rs[j] = log(At->a[At->cp[i]+j]+shift);
            ms[j] = At->ri[At->cp[i]+j];
        }
        for (j=0; j<nrats; j++) {
            mi=ms[j];
            ri=rs[j];
            for (k=0; k<nrats; k++) {
                if (k==j) { continue; } /* skip the same rating */
                mj=ms[k];
                rj=rs[k];
                Y[nitems*mj+mi] += ri-rj;
                N[nitems*mj+mi]++;
            }
        }
    }
    free(ms);
    free(rs);
}

void skew_binary(cs *At, double* Y, double* N, double equiv) {
    mwSize nusers = At->n;
    mwSize nitems = At->m;
    mwIndex i, j, k, mi, mj;
    double ri, rj;
    /* allocate memory for local caching, it seems to be faster */
    mwIndex *ms = calloc(nitems, sizeof(mwIndex));
    double *rs = calloc(nitems, sizeof(double));
    for (i=0; i<nusers; i++) {
        /* examine column i of At, which has the ratings for user i */
        mwIndex nrats = At->cp[i+1]-At->cp[i];
        for (j=0; j<nrats; j++) {
            rs[j] = At->a[At->cp[i]+j];
            ms[j] = At->ri[At->cp[i]+j];
        }
        for (j=0; j<nrats; j++) {
            mi=ms[j];
            ri=rs[j];
            for (k=0; k<nrats; k++) {
                if (k==j) { continue; } /* skip the same rating */
                rj=rs[k];
                mj=ms[k];
                if (ri>(rj+equiv)) {
                    Y[nitems*mj+mi] += 1.;
                } else if (rj>(ri+equiv)) {
                    Y[nitems*mj+mi] -= 1.;
                }
                N[nitems*mj+mi]++;
            }
        }
    }
    free(ms);
    free(rs);
}

void skew_binary_notie(cs *At, double* Y, double* N, double equiv) {
    mwSize nusers = At->n;
    mwSize nitems = At->m;
    mwIndex i, j, k, mi, mj;
    
    double ri, rj;
    /* allocate memory for local caching, it seems to be faster */
    mwIndex *ms = calloc(nitems, sizeof(mwIndex));
    double *rs = calloc(nitems, sizeof(double));
    for (i=0; i<nusers; i++) {
        /* examine column i of At, which has the ratings for user i */
        mwIndex nrats = At->cp[i+1]-At->cp[i];
        for (j=0; j<nrats; j++) {
            rs[j] = At->a[At->cp[i]+j];
            ms[j] = At->ri[At->cp[i]+j];
        }
        for (j=0; j<nrats; j++) {
            mi=ms[j];
            ri=rs[j];
            for (k=0; k<nrats; k++) {
                if (k==j) { continue; } /* skip the same rating */
                rj=rs[k];
                if (fabs(ri-rj)<=equiv) { continue; }
                mj=ms[k];
                Y[nitems*mj+mi] += 2.*(ri>rj)-1.;
                N[nitems*mj+mi]++;
            }
        }
    }
    free(ms);
    free(rs);    
}


void skew_logodds(cs *At, double* Y, double* N, double equiv, int pseudo) {
    mwSize nusers = At->n;
    mwSize nitems = At->m;
    mwIndex i, j, k, mi, mj;
    double ri, rj;
    /* allocate memory for local caching, it seems to be faster */
    mwIndex *ms = calloc(nitems, sizeof(mwIndex));
    double *rs = calloc(nitems, sizeof(double));
    for (i=0; i<nusers; i++) {
        /* examine column i of At, which has the ratings for user i */
        mwIndex nrats = At->cp[i+1]-At->cp[i];
        for (j=0; j<nrats; j++) {
            rs[j] = At->a[At->cp[i]+j];
            ms[j] = At->ri[At->cp[i]+j];
        }
        for (j=0; j<nrats; j++) {
            mi=ms[j];
            ri=rs[j];
            for (k=0; k<nrats; k++) {
                if (k==j) { continue; } /* skip the same rating */
                rj=rs[k];
                mj=ms[k];
                if (ri>=(rj+equiv)) {
                    Y[nitems*mj+mi] += 1.;
                } 
                N[nitems*mj+mi]++;
            }
        }
    }
    free(ms);
    free(rs);
    
    /* set lower triangular part to the correct values */
    for (j=0; j<nitems; j++) {
        for (i=j+1; i<nitems; i++) {
            if (N[nitems*j+i] > 0) {
                double Yji = Y[nitems*j+i];
                double Yij = Y[nitems*i+j];
                Y[nitems*j+i] = log(
                                    (Yji+(double)pseudo)/(Yij+(double)pseudo)
                                )*N[nitems*j+i];
                Y[nitems*i+j] = log(
                                    (Yij+(double)pseudo)/(Yji+(double)pseudo)
                                )*N[nitems*j+i];
            }
        }
    }
}

        
void mexFunction(
        int nargout, mxArray *pargout[],
        int nargin, const mxArray *pargin[])
{
    cs At;
    int type, pseudo;
    double shift, equiv=0.;
    double* Y, *N;
    if (nargin != 5 | nargout != 2) {
        mexErrMsgTxt(
        " Usage: [Y,N] = ratings2skew_mex(A,type,shift,equiv,pseudo); ");
    }
    cs_mex_get_sparse(pargin[0],&At);
    type = (int)mxGetScalar(pargin[1]);
    shift = mxGetScalar(pargin[2]);
    equiv = mxGetScalar(pargin[3]);
    pseudo = (int)mxGetScalar(pargin[4]);
    
    pargout[0] = mxCreateDoubleMatrix(At.m,At.m,mxREAL);
    Y = mxGetPr(pargout[0]);
    pargout[1] = mxCreateDoubleMatrix(At.m,At.m,mxREAL);
    N = mxGetPr(pargout[1]);    
    switch (type) {
        case 1: skew_arithmean(&At,Y,N); break;
        case 2: skew_geomean(&At,Y,N,shift); break;
        case 3: skew_binary(&At,Y,N,equiv); break;
        case 4: skew_binary_notie(&At,Y,N,equiv); break;
        case 5: skew_logodds(&At,Y,N,equiv,pseudo); break;
        default: mexErrMsgTxt("Invalid type!"); break;
    }
}

