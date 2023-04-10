/*************************************************************************
 * MATLAB MEX ROUTINE My3DMat2Vec.c
 *
 * [b3vec,xvec,yvec,zvec] = My3DMat2Vec(mat,xg,yg,zg)
 *
 ************************************************************************/

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    mwSize ndims, *dims;
    long int IsLast, totnele, dimz, dimx, dimy, dd,  ii, jj,  kk, cnt;
    double *mat,*xg,*yg,*zg,*b3vec,*xvec,*yvec,*zvec;
    
    mat=mxGetPr(prhs[0]);
    xg=mxGetPr(prhs[1]);
    yg=mxGetPr(prhs[2]);
    zg=mxGetPr(prhs[3]);
    
    ndims = mxGetNumberOfDimensions(prhs[0]);
    dims = mxGetDimensions(prhs[0]);
    
    dimz = (long int)dims[0];
    dimy = (long int)dims[1]; 
    dimx = (long int)dims[2]; 
    totnele=dimz*dimy*dimx;
    
    plhs[0]=mxCreateDoubleMatrix(totnele,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(totnele,1,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(totnele,1,mxREAL);
    plhs[3]=mxCreateDoubleMatrix(totnele,1,mxREAL);
    
    b3vec=mxGetPr(plhs[0]);
    xvec=mxGetPr(plhs[1]);
    yvec=mxGetPr(plhs[2]);
    zvec=mxGetPr(plhs[3]);
    
    cnt=0;
    for(ii=0; ii<dimx; ii=ii+1){
        for(jj=0; jj<dimy; jj=jj+1){
            for(kk=0; kk<dimz; kk=kk+1){                
                b3vec[cnt]=mat[ii*dimy*dimz+jj*dimz+kk];
                xvec[cnt]=xg[ii];
                yvec[cnt]=yg[jj];
                zvec[cnt]=zg[kk];
                cnt=cnt+1;
            }
        }
    }
} 