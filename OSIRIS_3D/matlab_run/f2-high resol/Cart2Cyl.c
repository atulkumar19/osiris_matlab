/*************************************************************************
 * MATLAB MEX ROUTINE Cart2Cyl.c
 *
 * [b1cyl,b2cyl,b3cyl]=Cart2Cyl(b1cart,b2cart,b3cart,x1cart,x2cart,x3cart)
 *
 ************************************************************************/

#include "mex.h"
#include "matrix.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    mwSize ndims, *dims;
    long int dimx, dimy, dimz,  ii, jj,  kk;
    double *b1cart,*b2cart,*b3cart,*x1cart,*x2cart,*x3cart;
    double *b1cyl,*b2cyl,*b3cyl;
    
    long int ind,iidimydimz,jjdimz;
    double xtmp,ytmp,rtmp,thetatmp,costheta,sintheta;
    
    b1cart=mxGetPr(prhs[0]);
    b2cart=mxGetPr(prhs[1]);
    b3cart=mxGetPr(prhs[2]);
    x1cart=mxGetPr(prhs[3]);
    x2cart=mxGetPr(prhs[4]);
    x3cart=mxGetPr(prhs[5]);
    
    ndims = mxGetNumberOfDimensions(prhs[0]);
    dims = mxGetDimensions(prhs[0]);
    
    dimz = (long int)dims[0];
    dimy = (long int)dims[1]; 
    dimx = (long int)dims[2]; 
        
    plhs[0]=mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[1]=mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[2]=mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
            
    b1cyl=mxGetPr(plhs[0]);
    b2cyl=mxGetPr(plhs[1]);
    b3cyl=mxGetPr(plhs[2]);
        
    for(ii=0; ii<dimx; ii=ii+1){
        iidimydimz=ii*dimy*dimz;
        xtmp=x1cart[ii];
        for(jj=0; jj<dimy; jj=jj+1){
            jjdimz=jj*dimz;
            ytmp=x2cart[jj];
            /*
            rtmp=sqrt(pow(xtmp,2)+pow(ytmp,2));
             */
            thetatmp=atan2(ytmp,xtmp);
            mexPrintf("theta=%lf\n",thetatmp);
            costheta=cos(thetatmp);
            sintheta=sin(thetatmp);
            for(kk=0; kk<dimz; kk=kk+1){
                ind=iidimydimz+jjdimz+kk;
                b1cyl[ind]=costheta*b1cart[ind]+sintheta*b2cart[ind];
                b2cyl[ind]=-sintheta*b1cart[ind]+costheta*b2cart[ind];
                b3cyl[ind]=b3cart[ind];
            }
        }
    }
} 