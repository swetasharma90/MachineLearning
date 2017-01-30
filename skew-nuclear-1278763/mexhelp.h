/**
 * Helper functions for mex files.  These functions automatic and
 * systematize operations that are common in mex files, such as checking
 * the input for the correct dimensions.
 * 
 * Many of the routines are taken from CSparse
 *
 * @author David F. Gleich or Tim Davis (taken from CSparse)
 *
 * History
 * -------
 * :2009-11-08: Initial coding
 */

#ifdef MATLAB_MEX_FILE
#define malloc mxMalloc
#define free mxFree
#define realloc mxRealloc
#define calloc mxCalloc
#endif

#define CS_MAX(a,b) ((a)>(b) ? (a) : (b))
/* wrapper for malloc */
void *cs_malloc (int n, size_t size)
{
    return (malloc (CS_MAX (n,1) * size)) ;
}

/* wrapper for calloc */
void *cs_calloc (int n, size_t size)
{
    return (calloc (CS_MAX (n,1), size)) ;
}

/* wrapper for free */
void *cs_free (void *p)
{
    if (p) free (p) ;       /* free p if it is not already NULL */
    return (NULL) ;         /* return NULL to simplify the use of cs_free */
}



/* check MATLAB input argument */
void cs_mex_check (int nel, int m, int n, int square, int sparse, int values,
    const mxArray *A)
{
    int nnel, mm = mxGetM (A), nn = mxGetN (A) ;
    if (values)
    {
        if (mxIsComplex (A))
        {
            mexErrMsgTxt ("matrix must be real; try CXSparse instead") ;
        }
    }
    if (sparse && !mxIsSparse (A)) mexErrMsgTxt ("matrix must be sparse") ;
    if (!sparse)
    {
        if (mxIsSparse (A)) mexErrMsgTxt ("matrix must be full") ;
        if (values && !mxIsDouble (A)) mexErrMsgTxt ("matrix must be double") ;
    }
    if (nel)
    {
        /* check number of elements */
        nnel = mxGetNumberOfElements (A) ;
        if (m >= 0 && n >= 0 && m*n != nnel) mexErrMsgTxt ("wrong length") ;
    }
    else
    {
        /* check row and/or column dimensions */
        if (m >= 0 && m != mm) mexErrMsgTxt ("wrong dimension") ;
        if (n >= 0 && n != nn) mexErrMsgTxt ("wrong dimension") ;
    }
    if (square && mm != nn) mexErrMsgTxt ("matrix must be square") ;
}

/* get a MATLAB int32 vector */
int *cs_mex_get_int_data(mwSize n, const mxArray* X) 
{
    if (n != mxGetNumberOfElements (X) ) {
        mexErrMsgTxt("wrong length");
    }
    if (mxGetClassID(X) != mxINT32_CLASS) {
        mexErrMsgTxt("matrix must be int32");
    }
    return (mxGetData(X)) ;
}

/* get a MATLAB uint32 vector */
unsigned int *cs_mex_get_uint_data(mwSize n, const mxArray* X) 
{
    if (n != mxGetNumberOfElements (X) ) {
        mexErrMsgTxt("wrong length");
    }
    if (mxGetClassID(X) != mxUINT32_CLASS) {
        mexErrMsgTxt("matrix must be uint32");
    }
    return (mxGetData(X)) ;
}


/* get a MATLAB dense column vector */
double *cs_mex_get_double (int n, const mxArray *X)
{
    cs_mex_check (1, n, 1, 0, 0, 1, X) ;
    return (mxGetPr (X)) ;
}

double *cs_mex_get_double_mat (int m, int n, const mxArray *X)
{
    cs_mex_check (0, m, n, 0, 0, 1, X) ;
    return (mxGetPr (X)) ;
}

/* return a double vector to MATLAB */
double *cs_mex_put_double (int n, const double *b, mxArray **X)
{
    double *x ;
    int k ;
    *X = mxCreateDoubleMatrix (n, 1, mxREAL) ;      /* create x */
    x = mxGetPr (*X) ;
    for (k = 0 ; k < n ; k++) x [k] = b [k] ;       /* copy x = b */
    return (x) ;
}

/* get a MATLAB flint array and convert to int */
int *cs_mex_get_int (int n, const mxArray *Imatlab, int *imax, int lo)
{
    double *p ;
    int i, k, *C = cs_malloc (n, sizeof (int)) ;
    cs_mex_check (1, n, 1, 0, 0, 1, Imatlab) ;
    p = mxGetPr (Imatlab) ;
    *imax = 0 ;
    for (k = 0 ; k < n ; k++)
    {
        i = p [k] ;
        C [k] = i - 1 ;
        if (i < lo) mexErrMsgTxt ("index out of bounds") ;
        *imax = CS_MAX (*imax, i) ;
    }
    return (C) ;
}

/* return an int array to MATLAB as a flint row vector */
mxArray *cs_mex_put_int (int *p, int n, int offset, int do_free)
{
    mxArray *X = mxCreateDoubleMatrix (1, n, mxREAL) ;
    double *x = mxGetPr (X) ;
    int k ;
    for (k = 0 ; k < n ; k++) x [k] = (p ? p [k] : k) + offset ;
    if (do_free) cs_free (p) ;
    return (X) ;
}


/* get a MATLAB flint array and convert to int */
mwIndex *cs_mex_get_index (mwSize n, const mxArray *Imatlab, 
            mwIndex *imax, mwIndex lo)
{
    double *p ;
    mwIndex i, k, *C = cs_malloc (n, sizeof (mwIndex)) ;
    cs_mex_check (1, n, 1, 0, 0, 1, Imatlab) ;
    p = mxGetPr (Imatlab) ;
    *imax = 0 ;
    for (k = 0 ; k < n ; k++)
    {
        i = p [k] ;
        C [k] = i - 1 ;
        if (i < lo) mexErrMsgTxt ("index out of bounds") ;
    }
    return (C) ;
}

/* return an int array to MATLAB as a flint row vector */
mxArray *cs_mex_put_index (mwIndex *p, mwSize n, mwIndex offset, 
            int do_free)
{
    mxArray *X = mxCreateDoubleMatrix (1, n, mxREAL) ;
    double *x = mxGetPr (X) ;
    mwIndex k ;
    for (k = 0 ; k < n ; k++) x [k] = (p ? p [k] : k) + offset ;
    if (do_free) cs_free (p) ;
    return (X) ;
}
