/*************************************************************************
 * MATLAB MEX ROUTINE fast_ind2sub.c
 *
 * [I1,I2,I3,...,In] = fast_ind2sub(SIZ,IND)
 * Similar to MATLAB ind2sub()
 * NOTES:
 * SIZ and IND must be double
 * returned results are uint32
 *
 * >> mex -O fast_ind2sub.c
 *
 * Author: Bruno Luong <brunoluong@xxxx.xxx>
 * History
 * Original: 22/April/2015
 ************************************************************************/

#include "mex.h"
#include "matrix.h"

#define SIZ prhs[0]
#define IND prhs[1]


/* Define correct type depending on platform */
#if defined(_MSC_VER) || defined(__BORLANDC__)
typedef unsigned __int32 uint32;
typedef signed __int32 int32;
typedef unsigned __int64 ulong64;
typedef signed __int64 long64;
#define INTMIN64 0x8000000000000000
#else
typedef unsigned int uint32;
typedef signed int int32;
typedef unsigned long long ulong64;
typedef signed long long long64;
#define INTMIN64 0x8000000000000000ULL
#endif

       
#define MAX_DIM 16

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    size_t nsubidx, n, m, i, j;
    mwSize ndims, *dims;
    mxArray *Temp;
    uint32 *ind, *ai, szi, indj;
    mxClassID ClassIND, ClassSIZ;
    double *ind_d, *size_d;
    int IsLast;
    
    nsubidx = nlhs;
    if (nsubidx < 1) nsubidx = 1;
    
    ndims = mxGetNumberOfDimensions(IND);
    dims = mxGetDimensions(IND);
    n = mxGetNumberOfElements(IND);
    m = mxGetNumberOfElements(SIZ);
    ClassIND = mxGetClassID(IND);
    ClassSIZ = mxGetClassID(SIZ);
    
    Temp = mxCreateNumericArray(ndims, dims, mxUINT32_CLASS, mxREAL);
    ind = (uint32*)mxGetData(Temp);

    if (ClassIND == mxDOUBLE_CLASS) {
        ind_d = (double*)mxGetData(IND);
    } else {
         mexErrMsgTxt("Class of input IND not supported");
    }

    if (ClassSIZ == mxDOUBLE_CLASS) {
        size_d = (double*)mxGetData(SIZ);
    } else {
         mexErrMsgTxt("Class of input SZ not supported");
    }
    
    for (j = 0; j < n; j++)
        ind[j] = (uint32)ind_d[j] - 1;
    
    for (i = 0; i < nsubidx; i++)
    {
        plhs[i] = mxCreateNumericArray(ndims, dims, mxUINT32_CLASS, mxREAL);
        ai = (uint32*)mxGetData(plhs[i]);
        IsLast = i == nsubidx-1;
        if (i < m)
        {
            szi = (uint32)size_d[i];
            if (IsLast)
                for (j = i+1; j < m; j++) szi *= (uint32)size_d[j];
        }
        else szi = 1;
        if (IsLast)
            /* optimization, update IND is not needed for the last subindice */
            for (j = 0; j < n; j++)
                ai[j] = 1 + (ind[j] % szi);
        else
            for (j = 0; j < n; j++)
            {
                indj = ind[j];
                ai[j] = 1 + (indj % szi);
                ind[j] = indj / szi;
            }
         
    }
    
    mxDestroyArray(Temp);
    
} 