/*=================================================================
 *
 * stackgdgradmex.c
 *
 * Requires a C99 compiler.
 *
 *=================================================================*/


#include <math.h>
#include <complex.h>
#include <string.h>
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

// Definitions
#define C 0.2997924580  // speed of light (um/fs)


// Core computation routine
static void stackgdgrad(
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
   double gd[],
   double r2grad[],
   double gdgrad[])
{
   // Variables
   int dx, kx;
   double neffs[3], dneffs[3], dterms[3], kneffs[3];
   double knefflay[n+1];
   complex dTdkgradterm, dTdkgradmat1[n+1], dTdkgradmat2[n+1];
   complex dTdkgradfor1, dTdkgradfor2;
   double pps[4], pms[4], dps[4];
   const int lastx = (n+1) & 1;  // last material index index
   int nx, px;  // material data indices
   complex ephi;
   complex Dlay[n+1], Tlay1[n+1], Tlay2[n+1], dTlay1[n+1], dTlay2[n+1];
   complex Tfor1[n+1], Tfor2[n+1], dTfor1[n+1], dTfor2[n+1];
   complex Trev1[n+1], Trev2[n+1], dTrev1[n+1], dTrev2[n+1];
   complex dTgradfor1[n], dTgradfor2[n];
   complex T1, T2;  // TODO: move?
   complex R, dR;
   double r, dr, dphi;
   complex Tgrad1, Tgrad2, dTgrad1, dTgrad2;
   double rgrad, phigrad, dphigrad;
   complex Rgrad, dRgrad;

   // Common terms and offset vectors.
   const double n0sintheta2 = (n0*sin(theta))*(n0*sin(theta));  // dispersionless
   const double *nsofk = ns, *dnsofk = dns;  // offset index pointers
   double *r2gradofk = r2grad, *gdgradofk = gdgrad;

   // Initialization.
   dTlay1[n] = 0; dTlay2[n] = 0;
   dTrev1[n] = 0; dTrev2[n] = 0;

   /*
    * Loop over all wavenumbers.
    */
   for (kx = 0; kx < nk; kx++)
   {
      /*
       * Precalculate all material parameters.
       * index variables: [n1, n2, nsub]
       * p variables: [p01, p12, p21, p2sub]
       */
      {
	 int j;
	 double pTEs[4], dpTEs[4];

	 // Calculate effective indices and common expressions.
	 for (j = 0; j < 3; j++)
	 {
	    neffs[j] = sqrt(nsofk[j]*nsofk[j] - n0sintheta2);
	    dneffs[j] = nsofk[j]*dnsofk[j]/neffs[j];
	    dterms[j] = neffs[j] + ks[kx]*dneffs[j];  // FIXME: neccesary?
	    kneffs[j] = ks[kx]*neffs[j];
	 }
	 // Calculate pTEs (needed regardless of polarization).
	 pTEs[0] = n0*cos(theta)/neffs[0];  // neff0/neff1
	 pTEs[1] = neffs[0]/neffs[1];
	 pTEs[2] = 1/pTEs[1];
	 pTEs[3] = neffs[lastx]/neffs[2];
	 dpTEs[0] = -pTEs[0]*dneffs[0]/neffs[0];
	 dpTEs[1] = (dneffs[0] - pTEs[1]*dneffs[1])/neffs[1];
	 dpTEs[2] = (dneffs[1] - pTEs[2]*dneffs[0])/neffs[0];
	 dpTEs[3] = (dneffs[lastx] - pTEs[3]*dneffs[2])/neffs[2];

	 // Calculate TM or TE as needed.
	 if (isTM)  // TM
	 {
	    double ps, p0s[4], dp0s[4];

	    // Calculate pTM from p0 (normal incidence) and pTE.
	    p0s[0] = n0/nsofk[0];
	    p0s[1] = nsofk[0]/nsofk[1];
	    p0s[2] = 1/p0s[1];
	    p0s[3] = nsofk[lastx]/nsofk[2];
	    dp0s[0] = -p0s[0]*dnsofk[0]/nsofk[0];
	    dp0s[1] = (dnsofk[0] - p0s[1]*dnsofk[1])/nsofk[1];
	    dp0s[2] = (dnsofk[1] - p0s[2]*dnsofk[0])/nsofk[0];
	    dp0s[3] = (dnsofk[lastx] - p0s[3]*dnsofk[2])/nsofk[2];
	    for (j = 0; j < 4; j++)
	    {
	       ps = pTEs[j]/(p0s[j]*p0s[j]);
	       dps[j] = (p0s[j]*dpTEs[j]-2*pTEs[j]*dp0s[j])/(p0s[j]*p0s[j]*p0s[j]);
	       pps[j] = 1 + ps;
	       pms[j] = 1 - ps;
	    }
	 }
	 else  // TE
	 {
	    for (j = 0; j < 4; j++)
	    {
	       dps[j] = dpTEs[j];
	       pps[j] = 1 + pTEs[j];
	       pms[j] = 1 - pTEs[j];
	    }
	 }
      }

      /*
       * Precalculate all layer matrices.
       */
      for (dx = 0; dx < n; dx++)
      {
	 // Select appropriate material parameter index
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
	 ephi = (cos(ds[dx]*kneffs[nx]) - I*sin(ds[dx]*kneffs[nx]))/2;
	 Dlay[dx] = -I*ds[dx]*dterms[nx];  // differential operator
	 Tlay1[dx] = ephi*pps[px];  // layer T element 1
	 Tlay2[dx] = ephi*pms[px];  // layer T element 2
	 dTlay1[dx] = Dlay[dx]*Tlay1[dx] + ephi*dps[px];
	 dTlay2[dx] = Dlay[dx]*Tlay2[dx] - ephi*dps[px];
	 knefflay[dx] = -kneffs[nx];
	 dTdkgradterm = -(ds[dx]*kneffs[nx] + I)*(neffs[nx] + ks[kx]*dneffs[nx]);
	 dTdkgradmat1[dx] = dTdkgradterm - I*kneffs[nx]*dps[px]/pps[px];
	 dTdkgradmat2[dx] = dTdkgradterm + I*kneffs[nx]*dps[px]/pms[px];
      }
      // Propagation into substrate matrix.
      Tlay1[n] = pps[3]/2;
      Tlay2[n] = pms[3]/2;
      dTlay1[n] = dps[3]/2;
      dTlay2[n] = -dps[3]/2;
      
      /*
       * Step forward through structure, calculating forward matrices.
       */
      Tfor1[0] = Tlay1[0]; Tfor2[0] = Tlay2[0];
      dTfor1[0] = dTlay1[0]; dTfor2[0] = dTlay2[0];
      dTgradfor1[0] = (1/ds[0] + I*knefflay[0])*dTfor1[0];
      dTgradfor2[0] = (1/ds[0] + I*knefflay[0])*dTfor2[0];
      for (dx = 1; dx < n+1; dx++)  // [2:n]
      {
	 // Calculate Lth forward matrix.
	 Tfor1[dx] = Tlay1[dx]*Tfor1[dx-1] + Tlay2[dx]*conj(Tfor2[dx-1]);
	 Tfor2[dx] = Tlay1[dx]*Tfor2[dx-1] + Tlay2[dx]*conj(Tfor1[dx-1]);

	 // Calculate total k derivative of the Lth forward T matrix.
	 dTfor1[dx] = dTlay1[dx]*Tfor1[dx-1] + dTlay2[dx]*conj(Tfor2[dx-1]) +
	    Tlay1[dx]*dTfor1[dx-1] + Tlay2[dx]*conj(dTfor2[dx-1]);
	 dTfor2[dx] = dTlay1[dx]*Tfor2[dx-1] + dTlay2[dx]*conj(Tfor1[dx-1]) +
	    Tlay1[dx]*dTfor2[dx-1] + Tlay2[dx]*conj(dTfor1[dx-1]);
      }

      /*
       * Step backward through structure, calculating reverse matrices.
       */
      Trev1[n] = Tlay1[n]; Trev2[n] = Tlay2[n];
      dTrev1[n] = dTlay1[n]; dTrev2[n] = dTlay2[n];
      for (dx = n-1; dx > 0; dx--)  // [n:-1:2]
      {
	 Trev1[dx] = Trev1[dx+1]*Tlay1[dx] + Trev2[dx+1]*conj(Tlay2[dx]);
	 Trev2[dx] = Trev1[dx+1]*Tlay2[dx] + Trev2[dx+1]*conj(Tlay1[dx]);

	 dTrev1[dx] = dTrev1[dx+1]*Tlay1[dx] + dTrev2[dx+1]*conj(Tlay2[dx]) +
	    Trev1[dx+1]*dTlay1[dx] + Trev2[dx+1]*conj(dTlay2[dx]);
	 dTrev2[dx] = dTrev1[dx+1]*Tlay2[dx] + dTrev2[dx+1]*conj(Tlay1[dx]) +
	    Trev1[dx+1]*dTlay2[dx] + Trev2[dx+1]*conj(dTlay1[dx]);
      }

      /*
       * Calculate full stack scalars at the current wavelength.
       */
      T1 = Tfor1[n];
      T2 = Tfor2[n];
      R = -T2/T1;
      r = cabs(R);
      r2[kx] = r*r;
      dR = (T2*dTfor1[n] - T1*dTfor2[n])/(T1*T1);
      dr = (creal(dR)*creal(R) + cimag(dR)*cimag(R))/r;
      dphi = (cimag(dR)*creal(R) - creal(dR)*cimag(R))/r2[kx];
      gd[kx] = dphi/C;

      /*
       * Step through structure, calculating gradients and output.
       */
      for (dx = 0; dx < n; dx++)  // [1:n-1]
      {
	 Tgrad1 = I*knefflay[dx]*
	    (Trev1[dx+1]*Tfor1[dx] - Trev2[dx+1]*conj(Tfor2[dx]));
	 Tgrad2 = I*knefflay[dx]*
	    (Trev1[dx+1]*Tfor2[dx] - Trev2[dx+1]*conj(Tfor1[dx]));
	 
	 // TODO: this could be done by having Tfor1[dx-1] = 1 and
	 // dTfor1[dx-1] = 0. That would eliminate the branch, which
	 // gets called nk*nl times even though it's triggered only
	 // once. The other option would be to try branch prediction hints.
	 if (dx > 0)
	 {
	    dTdkgradfor1 = dTdkgradmat1[dx]*Tlay1[dx]*Tfor1[dx-1] +
	       dTdkgradmat2[dx]*Tlay2[dx]*conj(Tfor2[dx-1]);
	    dTdkgradfor1 += I*knefflay[dx]*(Tlay1[dx]*dTfor1[dx-1] +
	    Tlay2[dx]*conj(dTfor2[dx-1]));
	    dTdkgradfor2 = dTdkgradmat1[dx]*Tlay1[dx]*Tfor2[dx-1] +
	       dTdkgradmat2[dx]*Tlay2[dx]*conj(Tfor1[dx-1]);
	    dTdkgradfor2 += I*knefflay[dx]*(Tlay1[dx]*dTfor2[dx-1] +
	    Tlay2[dx]*conj(dTfor1[dx-1]));
	 }
	 else
	 {
	    dTdkgradfor1 = dTdkgradmat1[dx]*Tlay1[dx];
	    dTdkgradfor2 = dTdkgradmat2[dx]*Tlay2[dx];
	 }
	 dTgrad1 = Trev1[dx+1]*dTdkgradfor1 + Trev2[dx+1]*conj(dTdkgradfor2);
	 dTgrad2 = Trev1[dx+1]*dTdkgradfor2 + Trev2[dx+1]*conj(dTdkgradfor1);
	 dTgrad1 += I*knefflay[dx]*(dTrev1[dx+1]*Tfor1[dx] -
				    dTrev2[dx+1]*conj(Tfor2[dx]));
	 dTgrad2 += I*knefflay[dx]*(dTrev1[dx+1]*Tfor2[dx] -
				    dTrev2[dx+1]*conj(Tfor1[dx]));

	 // Final gradients
	 Rgrad = -(R*Tgrad1 + Tgrad2)/T1;
	 rgrad = (creal(Rgrad)*creal(R) + cimag(Rgrad)*cimag(R))/r;
	 r2gradofk[dx] = 2*r*rgrad;
	 dRgrad = (Tgrad2*dTfor1[n] + Tgrad1*dTfor2[n] +
		   R*(2*Tgrad1*dTfor1[n] - T1*dTgrad1) - T1*dTgrad2)/(T1*T1);
	 phigrad = (cimag(Rgrad)*creal(R) - creal(Rgrad)*cimag(R))/r2[kx];
	 dphigrad = (cimag(dRgrad)*creal(R) - creal(dRgrad)*cimag(R))/r2[kx] -
	    (phigrad*dr + rgrad*dphi)/r;
	 gdgradofk[dx] = dphigrad/C;
      }

      // Update offset vectors for next wavenumber.
      nsofk += 3;
      dnsofk += 3;
      r2gradofk += n;
      gdgradofk += n;
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

#ifdef DEBUG
   // Check for proper number of arguments.
   if (nrhs < 7)
   {
      mexErrMsgTxt("stackgdfast2mex: 7 input arguments required.");
   }
   else if (nlhs != 4)
   {
      mexErrMsgTxt("stackgdfast2mex: 4 output arguments required.");
   }

   // Check for correct input dimensions.
   int mk = mxGetM(KS_IN);
   int m = mxGetN(DS_IN);
   if ((mk != 1) || (m != 1))
   {
      mexErrMsgTxt("StackGDFast2 called with incorrect dimensions.");
   }
#endif

   // Create matrices for return arguments
   R2_OUT = mxCreateDoubleMatrix(1, nk, mxREAL);
   GD_OUT = mxCreateDoubleMatrix(1, nk, mxREAL);
   R2GRAD_OUT = mxCreateDoubleMatrix(n, nk, mxREAL);
   GDGRAD_OUT = mxCreateDoubleMatrix(n, nk, mxREAL);

   // Assign pointers and values to the parameters
   double* r2 = mxGetPr(R2_OUT);
   double* gd = mxGetPr(GD_OUT);
   double* r2grad = mxGetPr(R2GRAD_OUT);
   double* gdgrad = mxGetPr(GDGRAD_OUT);

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
	 mexWarnMsgTxt("Negative layer thickness.\n");
	 break;
      }
   }
#endif

   // Do the actual computation
   stackgdgrad(nk, n, ks, ds, n0, ns, dns, theta, isTM,
		r2, gd, r2grad, gdgrad);

   return;
}
