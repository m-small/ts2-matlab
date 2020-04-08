/*=================================================================
 *
 * NEARESTK.C	.MEX file corresponding to NEAREST.M
 *	        returns the index vector of nearest neighbours to 
 *              each embedded point in x
 *
 * The calling syntax is:
 *
 *		[ind] = nearestk(X,tau,k)
 *
 *
 * This is a MEX-file for MATLAB.  
 *=================================================================*/
/* $Revision: 1.5 $ */
#include "math.h"
#include "mex.h"

void nearest(double	*x,
	     double	*ind,
	     int kneigh, /*number of near neighbours to find*/
	     int tau,
	     int m,
	     int n)
/* for each column of x find the index of the RMS closest other column, excluding 
   those closer than tau to that column*/
{
  int i,j,k,oj,oi=0;
  int *closest;
  double *bestdist;  
  double diff,dist,worstbestdist;
  
  closest = (int *) malloc(kneigh*sizeof(int));
  bestdist = (double *) malloc(kneigh*sizeof(double));
  
  for (i=0; i<n; i++) /*for each column of the matrix x*/{

    for (k=0; k<kneigh; k++) {        

    	*(closest+k)=0;
    	*(bestdist+k)=1e32;
    	}
    	worstbestdist=1e32;
    	     

    oj = 0;
    for (j=0; j<n; j++) /*compare to every other column*/{

      if (abs(i-j)>tau) /*the exclusion zone*/ {
	/*dist=rms(x[:,i]-x[:,j])*/
  	dist=0;

	for (k=0; k<m; k++) /* calculate rms */ {
	  diff= (*(x+oi+k)) - (*(x+oj+k));
	  dist += diff*diff;
	} /*of for k */

	if (dist<=worstbestdist) /*is this the best?*/ {
	  /*insert new closest, in the correct place amoung the other close ones */
	  /*... and update */
	  k=kneigh-1;
	  while (k>0 && *(bestdist+k-1)>dist) {
		*(bestdist+k)=*(bestdist+k-1);
		*(closest+k)=*(closest+k-1);
		k--;
		}
		*(bestdist+k)=dist;
		*(closest+k)=j;	
	} /*of while (k>0 && *(bestdist+k-1)>dist)*/
	worstbestdist=*(bestdist+kneigh-1);
      } /* of if (dist<=worstbestdist) */
      
      oj += m; /* =j*m */

    } /* of for j */

	for (k=0; k<kneigh; k++) {
	    *(ind+i*kneigh+k)=*(closest+k)+1;
	    }
    oi += m; /* =i*m */

  } /* of for i */
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    double *x,*ind;
    double dtau;
    int i,tau,mrows,ncols,mvect,nvect,kneigh; 
    bool needfree=0;
    
    /* Check for proper number of arguments */
    
    if (nrhs>3) { 
	mexErrMsgTxt("Too many input arguments."); 
    } else if (nrhs == 0) {
        mexErrMsgTxt("Insufficient input arguments.");
    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
    } 
    
    /*Assign a pointer to the input matrix*/
    x = mxGetPr(prhs[0]);

    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */     
    mrows = mxGetM(prhs[0]); 
    ncols = mxGetN(prhs[0]);
    

    /* Get the size of the "forbidden zone" --- the second input arg. */
    if (nrhs>=2){
      dtau=mxGetScalar(prhs[1]);
      if (floor(dtau)!=dtau) {
	    mexWarnMsgTxt("Second input should be an integer.");
      }
      tau=(int)floor(dtau);
    } else {
      tau=0;
    }

    /*number of nearest neighbours to find --- the third input arg */
    if (nrhs==3) {
	  kneigh=mxGetScalar(prhs[2]);
	  if (floor(kneigh)!=kneigh) {
		    mexWarnMsgTxt("Second input should be an integer.");
      }
      kneigh=(int)floor(kneigh);
      } else {
      kneigh=1;
      }
  

    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(kneigh, ncols, mxREAL); 


    /* Assign pointer to the ouput matrix */ 
    ind = mxGetPr(plhs[0]);

        
    /* Do the actual computations in a subroutine */
    nearest(x,ind,kneigh,tau,mrows,ncols); 

    return;
    
}


