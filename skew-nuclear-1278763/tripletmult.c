/*
 * Compute a matrix vector product given a triplet form (i,j,v) of a matrix.
 *
 * @author David F. Gleich
 *
 * History
 * -------
 * :2009-11-21: Initial coding
 */

#include <mex.h>
#include "mexhelp.h"

void mexFunction(
        int nargout, mxArray *pargout[],
        int nargin, const mxArray *pargin[])
{
    mwIndex i;
    double *v, *x, *y;
    
    
    mwSize nzs;
    
    mwSize m, n, k, N;
    if (nargin  != 6 || nargout != 1) {
        mexErrMsgTxt ("Usage: y = tripletmult(i,j,v,x,m,n)") ;
    }
    nzs = mxGetNumberOfElements(pargin[0]);
    v = cs_mex_get_double(nzs, pargin[2]);
    m = (mwSize)mxGetScalar(pargin[4]);
    n = (mwSize)mxGetScalar(pargin[5]);
    x = cs_mex_get_double(n, pargin[3]);
    pargout[0] = mxCreateDoubleMatrix(m,1,mxREAL);
    y = mxGetPr(pargout[0]);
    
    if (mxGetClassID(pargin[0]) == mxINT32_CLASS) {
        int *ii=cs_mex_get_int_data(nzs,pargin[0]);
        int *ij=cs_mex_get_int_data(nzs,pargin[1]);
        for (i=0; i<nzs; i++) {
            y[ii[i]] += v[i]*x[ij[i]];
        }
    } else if (mxGetClassID(pargin[0]) == mxUINT32_CLASS) {
        unsigned int *ii=cs_mex_get_uint_data(nzs,pargin[0]);
        unsigned int *ij=cs_mex_get_uint_data(nzs,pargin[1]);
        for (i=0; i<nzs; i++) {
            y[ii[i]] += v[i]*x[ij[i]];
        }
    } else {
        double *fi, *fj;
        fi = cs_mex_get_double(nzs, pargin[0]);
        fj = cs_mex_get_double(nzs, pargin[1]);
        for (i=0; i<nzs; i++) {
            y[(mwIndex)(fi[i])-1] += v[i]*x[(mwIndex)(fj[i])-1];
        }
    }
}

