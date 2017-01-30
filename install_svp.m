%Script to compile mex files stored in "private" folder
%Written by: Prateek Jain (pjain@cs.utexas.edu) and Raghu Meka (raghu@cs.utexas.edu)
%Last updated: October, 2009

try
    cd private;
%      mex -O -v -lmwlapack -lmwblas bdsqr_mex.c dbdqr.f -UDBDQR_IN_C -output bdsqr
%     mex -O -v -lmwlapack -lmwblas reorth_mex.c reorth.f -output reorth
%     mex -O -v -lmwlapack -lmwblas tqlb_mex.c tqlb.f -output tqlb
    mex compute_X_Omega.c
    mex compute_Xss_Omega.c
    mex compOmegaYx.c
    mex compOmegaYtx.c
    cd ..
catch me
    cd ..
    rethrow(me)
end