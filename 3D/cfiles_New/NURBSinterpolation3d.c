#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "NURBS.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*	Return the NURBS basis function to matlab
	
	//
	// We expect the function to be called as 
	// [interp] = NURBSinterpolation(xi, eta, p, q, uknot, vknot
                                   points, weights)
	//
	//	xi   = point where we want to interpolate [xi eta]
	//	uknot = knot vector u direction
    //	vknot = knot vector v direction
	//	vector of points in format [pt1 pt2 pt3 pt4 .... ptn] */
	
    // Wang Zhetong [2024/5/29]

    // Validate the number of inputs
    if (nrhs != 8) {
        mexErrMsgTxt("You must pass in 8 arguments to the function. Usage: [interp]=NURBSinterpolation(xi, eta, p, q, uknot, vknot, points, weights)\n");
    }
    
    // Get the inputs
    double *xi      = mxGetPr(prhs[0]);    
    double *eta     = mxGetPr(prhs[1]);    
    double *p_in    = mxGetPr(prhs[2]);
    double *q_in    = mxGetPr(prhs[3]);
    double *uknot   = mxGetPr(prhs[4]);
    double *vknot   = mxGetPr(prhs[5]);
    double *points  = mxGetPr(prhs[6]);
    double *weight  = mxGetPr(prhs[7]);
    
    // Create output
    int dim = mxGetN(prhs[6]); // Number of columns in points (2 for 2D, 3 for 3D)
    if (dim != 2 && dim != 3) {
        mexErrMsgTxt("Points array must have 2 or 3 columns for 2D or 3D points.");
    }
    
    plhs[0] = mxCreateDoubleMatrix(1, dim, mxREAL); 
    double *interp  = mxGetPr(plhs[0]);
    
    // Convert inputs to appropriate types
    int p = (int) *p_in;
    int q = (int) *q_in;
    
    int numKnotU = mxGetN(prhs[4]);
    int numKnotV = mxGetN(prhs[5]);
    
    int mu = numKnotU - 1;
    int mv = numKnotV - 1;
    
    int nu = mu - p - 1;
    int nv = mv - q - 1;
    
    int numPoints = mxGetM(prhs[6]);
    int numWeights = mxGetN(prhs[7]);
    
    if (numPoints != numWeights) {
        mexErrMsgTxt("The number of points and weights must be the same.");
    }
    
    double tol = 100 * DBL_EPSILON;
    
    if (fabs(*xi - uknot[numKnotU - 1]) < tol) *xi  = uknot[numKnotU - 1] - tol;
    if (fabs(*eta - vknot[numKnotV - 1]) < tol) *eta = vknot[numKnotV - 1] - tol;
    
    // Allocate memory for basis functions
    double *N = (double *)malloc(sizeof(double) * (p + 1));
    double *M = (double *)malloc(sizeof(double) * (q + 1));
    
    int spanU = FindSpan(nu, p, *xi, uknot); 
    int spanV = FindSpan(nv, q, *eta, vknot); 
    
    BasisFuns(spanU, *xi, p, uknot, N);
    BasisFuns(spanV, *eta, q, vknot, M);
    
    // Initialize interpolation variables
    for (int i = 0; i < dim; i++) {
        interp[i] = 0.0;
    }
    double wght = 0.0;

    // Compute the approximation
    for (int l = 0; l <= q; l++) {
        double tempu = 0.0, tempv = 0.0, tempz = 0.0, tempw = 0.0;
        int vind = spanV - q + l;
        
        for (int k = 0; k <= p; k++) {
            int uind = spanU - p + k;
            int id = uind + vind * (nu+1); 
            
            tempu += N[k] * points[id] * weight[id];
            tempv += N[k] * points[id + numPoints] * weight[id];
            if (dim == 3) {
                tempz += N[k] * points[id + 2 * numPoints] * weight[id];
            }
            tempw += N[k] * weight[id];
        }
        
        interp[0] += tempu * M[l];
        interp[1] += tempv * M[l];
        if (dim == 3) {
            interp[2] += tempz * M[l];
        }
        wght += tempw * M[l];
    }
    
    // Projection
    for (int i = 0; i < dim; i++) {
        interp[i] /= wght;
    }
    
    // Free allocated memory
    free(N);
    free(M);
}
