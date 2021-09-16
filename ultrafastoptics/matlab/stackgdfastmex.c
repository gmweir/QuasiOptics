/*============================================================
 *
 * stackgdfastmex.c
 *
 * Computes gradients of reflectivity and group delay using
 * approximate analytic methods. The algorithm in this file is
 * implemented in a much more readable fashion in the corresponding M
 * file.
 * 
 * TODO: Try unions to save memory.
 * TODO: Make it work for odd numbers of layers! Idiot.
 * TODO: Allow for dispersive air.
 *
 *============================================================*/

// Includes
#include <math.h>
#include <complex.h>
#include <string.h>
#include <time.h>
#include "mex.h"

// Input Arguments
#define	KS_IN		prhs[0]
#define	DS_IN		prhs[1]
#define N0_IN		prhs[2]
#define NS_IN		prhs[3]
#define DNS_IN		prhs[4]
#define THETA_IN	prhs[5]
#define POL_IN		prhs[6]

// Output Arguments
#define	R2_OUT		plhs[0]
#define GD_OUT	        plhs[1]
#define R2GRAD_OUT      plhs[2]
#define GDGRAD_OUT      plhs[3]


// Definitions and Utility Functions
#define C 0.2997924580  // speed of light (um/fs)

struct tmatrix {
      complex a, b;
};

/* inline void tmatinit(struct tmatrix* y, complex v) */
/* { */
/*    y->a = v; */
/*    y->b = v; */
/* } */

#define TMATMUL(Y,X1,X2) Y.a = X1.a*X2.a + X1.b*conj(X2.b); \
Y.b = X1.a*X2.b + X1.b*conj(X2.a);
#define TMATSCAL(Y,S,X) Y.a = (S)*X.a; Y.b = (S)*X.b;
#define TMATSUM(Y,X1,X2) Y.a = X1.a + X2.a; Y.b = X1.b + X2.b;
#define TMATCOPY(Y,X) Y.a = X.a; Y.b = X.b;
#define TMATREVCOPY(Y,X) Y.a = X.a; Y.b = conj(X.b);
#define TMATREVMUL(Y,X1,X2) Y.a = X1.a*X2.a + conj(X1.b*X2.b); \
Y.b = X1.a*X2.b + conj(X1.b*X2.a);
#define TMATSCALREV(Y,S,X) Y.a = (S)*X.a; Y.b = (S)*conj(X.b);

//#define CREATEWATCH clock_t ticks1 = clock(); clock_t ticks2;
//#define HITWATCH(SECT) ticks2 = clock(); \
//mexPrintf("sec %d: %ld ticks.\n", SECT, ticks2 - ticks1); ticks1 = ticks2;


// *** Core computation routine ***
static void stackgdfast(
   const int nk,
   const int n,
   const double ks[],
   const double ds[],
   const double n0,
   const double ns[],
   const double dns[],
   const double theta,
   const int isTM,

   double r2[],
   double gd[])
{
   // Variables
   register int dx;  // layer array index
   double knefflay[n+1];  // common term
   complex Dlay[n+1], Play[n+1];  // derivative and propagation diags.
   struct tmatrix Rlay[n+1];  // reflection transfer matrix
   struct tmatrix dTlayTfor, TlaydTfor;
   struct tmatrix Tfor[n+1],  dTfor[n+1];
   struct tmatrix T;  // temp usage
   complex R, dR;
   double r, dr, dphi;

   // Common terms and offset vectors.
   const double n0sintheta2 = (n0*sin(theta))*(n0*sin(theta));
   const double *nsofk = ns, *dnsofk = dns;  // ptrs. to mat. indices

   /*
    * Loop over all wavenumbers.
    */
   for (int kx = 0; kx < nk; kx++)
   {
      /*
       * Precalculate all material parameters.
       * index array order: [n1, n2, nsub]
       * p array order: [p01, p12, p21, p2/1sub]
       */
      {
	 double neffs[3], dneffs[3];  // indices
	 int j;  // index for loops which should hopefully be unrolled
	 double dterms[3], kneffs[3];  // commonly used terms
	 double ps, pps[4], pms[4];  // index ratio factors
	 double pTEs[4];
	 int klast = (n+1)%2;  // index of last material (0/1)

	 // Calculate effective indices and common expressions.
	 for (j = 0; j < 3; j++)
	 {
	    neffs[j] = sqrt(nsofk[j]*nsofk[j] - n0sintheta2);
	    dneffs[j] = nsofk[j]*dnsofk[j]/neffs[j];
	    dterms[j] = neffs[j] + ks[kx]*dneffs[j];
	    kneffs[j] = ks[kx]*neffs[j];
	 }
	 // Calculate pTEs (needed regardless of polarization).
	 pTEs[0] = n0*cos(theta)/neffs[0];  // neff0/neff1
	 pTEs[1] = neffs[0]/neffs[1];
	 pTEs[2] = 1/pTEs[1];
	 pTEs[3] = neffs[klast]/neffs[2];
	 // Calculate TM or TE as needed.
	 if (isTM)  // TM
	 {
	    double p0s[4];

	    // Calculate pTM from p0 (normal incidence) and pTE.
	    p0s[0] = n0/nsofk[0];
	    p0s[1] = nsofk[0]/nsofk[1];
	    p0s[2] = 1/p0s[1];
	    p0s[3] = nsofk[klast]/nsofk[2];
	    for (j = 0; j < 4; j++)
	    {
	       ps = pTEs[j]/(p0s[j]*p0s[j]);
	       pps[j] = (1 + ps)/2;
	       pms[j] = (1 - ps)/2;
	    }
	 }
	 else  // TE
	 {
	    for (j = 0; j < 4; j++)
	    {
	       pps[j] = (1 + pTEs[j])/2;
	       pms[j] = (1 - pTEs[j])/2;
	    }
	 }

	 /*
	  * Precalculate all layer matrices.
	  */
	 int nx, px;  // material data array indices
	 for (dx = 0; dx < n; dx++)
	 {
	    // Select appropriate material parameters (by array index)
	    if (dx & 1) {  // even layer (mod(dx,2) == 1)
	       nx = 1;
	       px = 1; }
	    else if (dx == 0) {  // first layer
	       nx = 0;
	       px = 0; }
	    else {  // odd layer
	       nx = 0;
	       px = 2;
	    }

	    // ephi = exp(I*d*k*neff)/2
	    Dlay[dx] = -I*ds[dx]*dterms[nx];  // differential operator
	    Rlay[dx].a = pps[px];
	    Rlay[dx].b = pms[px];
	    Play[dx] = (cos(ds[dx]*kneffs[nx]) - I*sin(ds[dx]*kneffs[nx]));
	    knefflay[dx] = -kneffs[nx];  // $
	 }
	 // Propagation into substrate matrix.
	 Rlay[n].a = pps[3]; Rlay[n].b = pms[3];
      }  // material and layers block

      /*
       * Step forward through structure, calculating forward matrices.
       */
      TMATSCAL(Tfor[0], Play[0], Rlay[0]);
      TMATSCAL(dTfor[0], Dlay[0], Tfor[0]);
      for (dx = 1; dx < n; dx++)
      {
	 // Inductive T-matrix computations.
	 TMATSCAL(T, Play[dx], Rlay[dx]);
	 TMATMUL(Tfor[dx], T, Tfor[dx-1]);
	 TMATMUL(TlaydTfor, T, dTfor[dx-1]);
	 TMATSCAL(dTlayTfor, Dlay[dx], Tfor[dx]);
	 TMATSUM(dTfor[dx], dTlayTfor, TlaydTfor);
      }
      // Handle substrate.
      TMATMUL(Tfor[n], Rlay[n], Tfor[n-1]);
      TMATMUL(dTfor[n], Rlay[n], dTfor[n-1]);

      /*
       * Calculate full stack scalars at the current wavelength.
       */
      {
	 struct tmatrix T;

	 TMATCOPY(T, Tfor[n]);
	 R = -T.b/T.a;  // complex reflection
	 r = cabs(R);  // magnitude
	 r2[kx] = r*r;  // reflectivity
	 dR = (T.b*dTfor[n].a - T.a*dTfor[n].b)/(T.a*T.a);
	 dr = (creal(dR)*creal(R) + cimag(dR)*cimag(R))/r;
	 dphi = (cimag(dR)*creal(R) - creal(dR)*cimag(R))/r2[kx];
	 gd[kx] = dphi/C;
      }  // block

      // Update pointers for next wavenumber.
      nsofk += 3;
      dnsofk += 3;
   }

   return;
}


// MEX function gateway routine.
void mexFunction(
   int nlhs, mxArray* plhs[],
   int nrhs, const mxArray* prhs[])
{
   // Check dimensions.
   int nk = mxGetN(KS_IN);
   int n = mxGetM(DS_IN);

   // Check for proper number of arguments.
   if (nrhs < 7)
   {
      mexErrMsgTxt("stackgdfastmex: 7 input arguments required.");
   }
   else if (nlhs != 2)
   {
      mexErrMsgTxt("stackgdfastmex: 2 output arguments required.");
   }

   // Check for correct input dimensions.
   int mk = mxGetM(KS_IN);
   int m = mxGetN(DS_IN);
   if ((mk != 1) || (m != 1))
   {
      mexErrMsgTxt("stackgdfastmex called with incorrect dimensions.");
   }

   // Create matrices for return arguments
   R2_OUT = mxCreateDoubleMatrix(1, nk, mxREAL);
   GD_OUT = mxCreateDoubleMatrix(1, nk, mxREAL);

   // Assign pointers and values to the parameters
   double* r2 = mxGetPr(R2_OUT);
   double* gd = mxGetPr(GD_OUT);

   double* ks = mxGetPr(KS_IN);
   double* ds = mxGetPr(DS_IN);
   double n0 = mxGetScalar(N0_IN);
   double* ns = mxGetPr(NS_IN);
   double* dns = mxGetPr(DNS_IN);
   double theta = mxGetScalar(THETA_IN);
   int polstrlen = mxGetN(POL_IN);
   char polstr[8];
   mxGetString(POL_IN, polstr, polstrlen + 1);
   int isTM = (strncmp("TM", polstr, 2) == 0);

#ifdef DEBUG
   // Check for negative thicknesses.
   for (int dx = 0; dx < n; dx++)
   {
      if (ds[dx] < 0.0)
      {
	 mexWarnMsgTxt("stackgdfast: negative layer thickness.\n");
	 break;
      }
   }
#endif

   // Do the actual computation
   stackgdfast(nk, n, ks, ds, n0, ns, dns, theta, isTM,
	       r2, gd);

   return;
}
