#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NURBS.h"
#include <float.h>
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*	Return the NURBS basis function to matlab
	
	//
	// We expect the function to be called as [Nip] = NURBSBasis(i, p, xi, knot, weights)
	//
	// p = order of basis (0,1,2 etc)
	// knot = knot vector
	// i = basis function we want (1,2,3 ie. Hughes' notation)
	// xi = the coordinate we wish to evaluate at
	// weights = the weightings for the NURNS function (length(weights) = m - p -1) */
	
    // This version supports multiple basis function indices `i` as input.
    // It returns arrays of NURBS basis values and their derivatives at the specified parametric location `xi`.

	if(nrhs != 5) mexErrMsgTxt("You fool! You haven't passed in 5 arguments to the function."
								"We expect it to be in the form [Nip] = NURBSBasis(i, p, xi, knot, weights)\n");
								
	int c, k, p, m, i, numKnot, numWeights;
	
    int ni; // num of indices

	double *p_in, *m_in ,*i_in;
	double *knot, *xi, *weight;
	
	double tol=100*DBL_EPSILON;
	
	for(c = 0; c < nrhs; c++)
	{
		switch(c)
		{
			case 0:
				i_in = mxGetPr(prhs[c]);
                ni = mxGetN(prhs[c]);

			case 1:
				p_in = mxGetPr(prhs[c]);
				p = (int)*p_in;
				/*mexPrintf("\n\np=%d\n", p);*/

			case 2:
				xi = mxGetPr(prhs[c]);
				/*mexPrintf("\nxi=%2.20f\n", *xi);											*/

			case 3:
				knot = mxGetPr(prhs[c]);
				numKnot = mxGetN(prhs[c]);
				m = numKnot - 1;
				/*mexPrintf("\nWe have a knot vector with %d components\n", numKnot);
				for(k = 0; k < numKnot; k++) mexPrintf("%2.2f\t", *(knot+k));
				mexPrintf("\n");*/

			case 4:
				weight = mxGetPr(prhs[c]);
				numWeights = mxGetM(prhs[c]);

		}
	}	
	
	if(fabs(*xi-knot[numKnot-1]) < tol) 
		*xi = knot[numKnot-1] - tol;
	
	
	/* and call the basis function routine*/

	double w_interp, dw_interp_dxi ;
	
	int n = m - p -1;
	
	double *N = (double *)malloc(sizeof(double)*(p+1));
	double *dN = (double *)malloc(sizeof(double)*(p+1));
	double **ders = init2DArray(p+1, p+1);
	
	int span = FindSpan(n, p, *xi, knot);

	BasisFuns(span, *xi, p, knot, N);
	dersBasisFuns(span, *xi, p, p, knot, ders);	

	w_interp = 0.0; dw_interp_dxi = 0.0;
	
	for(c = 0; c <= p; c++)
	{
		w_interp += N[c] * weight[span-p+c];
		dw_interp_dxi += ders[1][c] * weight[span-p+c];
	}

    plhs[0] = mxCreateDoubleMatrix(1,ni,mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(1,ni,mxREAL); 
    
    double *Rip      = mxGetPr(plhs[0]);
    double *dRip_dxi  = mxGetPr(plhs[1]);
    double Nip;
    int temp;
    
    for (c = 0; c < ni; c++)
    {
        temp = (int)i_in[c] - 1;
        dersOneBasisFuns(p, m, knot, temp, *xi, p, dN);
        Nip = OneBasisFun(p, m, knot, temp, *xi);			
        
        Rip[c] = Nip * weight[temp] / w_interp;
		dRip_dxi[c] = weight[temp] * ( w_interp * dN[1] - dw_interp_dxi * Nip ) / (w_interp * w_interp);
    }

	free2Darray(ders, (p+1));
	free(N); free(dN);
}



 
 

