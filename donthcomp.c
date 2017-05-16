/* donthcomp.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

// #include "f2c.h"

// void c_donthcomp(double *ear, int ne, double * param, double *photar);

/* Table of constant values */

#include "common.h"

static double c_b2 = 2.;
static double c_b8 = 10.;
static double c_b10 = 4.;
static double c_b11 = 3.;
static double c_b42 = 1.;
static double c_b45 = .001;
static double c_b46 = -.66666666666666663;
static double c_b47 = 1.663753;
static double pi = 3.14159265358979323846;


/************ functions needed due to the Fortran-to-C conversion *****/

double pow_dd(double *ap, double *bp){
return(pow(*ap, *bp) );
}

static int int_min(int a, int b){
	return (a<b) ? a : b;
}

#define log10e 0.43429448190325182765
double d_lg10(double *x){
	return( log10e * log(*x) );
}

/*********************************************************************/


static int f_mcdint__(double *et, double *value)
{
    /* Initialized data */

    static double gc[3] = { .078196667,-1.066202,1.192418 };
    static double gw[3] = { .5207874,.513457,.4077983 };
    static double gn[3] = { .3728691,.039775528,.037766505 };
    static double res[98] = { 9.6198382e-4,.0010901181,.0012310012,
	    .0013841352,.0015481583,.0017210036,.0018988943,.002076939,
	    .0022484281,.0024049483,.0025366202,.0026316255,.0026774985,
	    .0026613059,.0025708784,.0023962965,.002130655,.0017725174,
	    .0013268656,8.0657672e-4,2.3337584e-4,-3.6291778e-4,-9.4443569e-4,
	    -.0014678875,-.0018873741,-.0021588493,-.0022448371,-.0021198179,
	    -.0017754602,-.0012246034,-5.0414167e-4,3.2507078e-4,.0011811065,
	    .0019673402,.0025827094,.0029342526,.0029517083,.0026012166,
	    .0018959062,9.0128649e-4,-2.6757144e-4,-.0014567885,-.002492855,
	    -.0032079776,-.0034678637,-.0031988217,-.0024080969,-.001193624,
	    2.6134145e-4,.0017117758,.0028906898,.0035614435,.0035711778,
	    .0028921374,.0016385898,4.9857464e-5,-.0015572671,-.0028578151,
	    -.0035924212,-.0036253044,-.002975086,-.0018044436,-3.7796664e-4,
	    .0010076215,.0020937327,.0027090854,.0028031667,.0024276576,
	    .0017175597,8.1030795e-4,-1.2592304e-4,-9.4888491e-4,-.0015544816,
	    -.0018831972,-.0019203142,-.0016905849,-.0012487737,-6.6789911e-4,
	    -2.7079461e-5,5.9931935e-4,.0011499748,.0015816521,.0018709224,
	    .0020129966,.0020184702,.0019089181,.0017122289,.001458377,
	    .0011760717,8.9046768e-4,6.2190822e-4,3.8553762e-4,1.9155022e-4,
	    4.5837109e-5,-4.9177834e-5,-9.3670762e-5,-8.9622968e-5,
	    -4.01538532e-5 };

    /* System generated locals */
    double d__1;


    /* Local variables */
    static int j;
    static double z__, pos, loget, gaufact, resfact;

    loget = d_lg10(et);
    pos = (loget - d_lg10(&c_b45)) / .06 + 1;
    j = (int) pos;
    if (j < 1) {
	resfact = res[0];
    } else if (j >= 98) {
	resfact = res[97];
    } else {
	pos -= j;
	resfact = res[j - 1] * (1. - pos) + res[j] * pos;
    }
    gaufact = 1.;
    for (j = 1; j <= 3; ++j) {
	z__ = (loget - gc[j - 1]) / gw[j - 1];
	gaufact += gn[j - 1] * exp(-z__ * z__ / 2.);
    }
    d__1 = *et / .001;
    *value = pow_dd(&d__1, &c_b46) * 193.21556 * (pow_dd(et, &c_b47) *
	    .52876731 + 1.) * exp(-(*et)) * gaufact * (resfact + 1.);
    return 0;
} /* f_mcdint__ */

/*     Multi-Color Disk SPECTRUM */
/*     SEE MITSUDA ET AL. 1984 PASJ 36, 741. */
/*     & MAKISHIMA ET AL. APJ 308, 635 1986) */
/*     NORMALIZATION={RIN(KM)/(D/10KPC)}^2*COS(THETA) */
/*     TIN=     INNER TEMPERATURE OF THE DISK */
/*     E  =     ENERGY */
/*     Rin2 = Normalization factor in the unit of above normalization. */
/*     PHOTON = PHOTONS/S/KEV/CM2 */
static int f_mcdspc__(double *e, double *tin, double *
	rin2, double *flux)
{
    static double et, value;

/*  E = X-ray energy (keV) */
/*  Tin = inner edge color-temperature (keV) */
/*  Rin2 = inner edge radius(km) ^2 * cos(inclination)/ [D (10kpc)]^2 */
/*  Flux = photon flux, photons/sec/cm^2/keV */
    if (*tin == 0.) {
	*flux = 0.;
	return 0;
    }
    et = *e / *tin;
    f_mcdint__(&et, &value);
    *flux = value * *tin * *tin * *rin2 / 361.;
    return 0;
} /* f_mcdspc__ */


static int f_xsdskb__(double *ear, int *ne, double *
	param, int *idt, double *photar, double *photer)
{
    /* Initialized data */

    static double gauss[10]	/* was [5][2] */ = { .236926885,.47862867,
	    .568888888,.47862867,.236926885,-.906179846,-.53846931,0.,
	    .53846931,.906179846 };

    /* System generated locals */
    int i__1;

    /* Local variables */
    static double e;
    static int i__, j;
    static double xh, xn, tin, photon;

/*     Multicolour disk blackbody model used in ISAS, Japan. */
/*     See Mitsuda et al. PASJ, 36, 741 (1984) */
/*     & Makishima et al. ApJ 308, 635 1986) */
/*     Ken Ebisawa 1992/12/22 */
/*     Modified to use double precision and to make numerical */
/*     integration faster. */
/*     Ken Ebisawa 1993/07/29 */
/*     Modified the algorithm. by Kazuhisa MITSUDA May 31, 1994 */
/*     A numerical calculation was first done within an accruacy of 0.01 %. */
/*     Then an interpolation formula was found. */
/*     The interpolation is precise within 1e-5 level. */
/*     These coefficients are taken from the spectral fitting program SPFD */
/*     in ISAS. */
/*     These are used for Gaussian Integral in the given energy band */
    /* Parameter adjustments */
    --photer;
    --photar;
    --param;

    /* Function Body */
/* suppress a warning message from the compiler */
    i__ = *idt;
/* this model has no errors */
    i__1 = *ne;
    for (i__ = 1; i__ <= i__1; ++i__) {
	photer[i__] = 0.f;
    }
    tin = param[1];
    i__1 = *ne;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xn = (ear[i__] - ear[i__ - 1]) / 2.f;
	photar[i__] = 0.f;
	xh = xn + ear[i__ - 1];
	for (j = 1; j <= 5; ++j) {
	    e = xn * gauss[j + 4] + xh;
	    f_mcdspc__(&e, &tin, &c_b42, &photon);
	    photar[i__] += (double) (gauss[j - 1] * photon);
/* L200: */
	}
	photar[i__] *= xn;
/* L100: */
    }
    return 0;
} /* f_xsdskb__ */


static int f_thermlc__(double *tautom, double *theta,
	double *deltal, double *x, int *jmax, double *dphesc,
	double *dphdot, double *bet, double *c2)
{
    /* System generated locals */
    int i__1;

    /* Builtin functions */
    double sqrt(double), pow_dd(double *, double *);

    /* Local variables */
    static double a[900], b[900], c__[900], d__[900], g[900];
    static int j;
    static double u[900], t1, t2, t3, w1, w2, aa, c20;
    static int jj;
    static double x32, gam[900], alp[900];

/* This program computes the effects of Comptonization by */
/* nonrelativistic thermal electrons in a sphere including escape, and */
/* relativistic corrections up to photon energies of 1 MeV. */
/* the dimensionless photon energy is x=hv/(m*c*c) */

/* The input parameters and functions are: */
/* dphdot(x), the photon production rate */
/* tautom, the Thomson scattering depth */
/* theta, the temperature in units of m*c*c */
/* c2(x), and bet(x), the coefficients in the K-equation and the */
/*   probability of photon escape per Thomson time, respectively, */
/*   including Klein-Nishina corrections */
/* The output parameters and functions are: */
/* dphesc(x), the escaping photon density */
/* u(x) is the dimensionless photon occupation number */
    /* Parameter adjustments */
    --c2;
    --bet;
    --dphdot;
    --dphesc;
    --x;

    /* Function Body */
    c20 = *tautom / *deltal;

/* determine u */
/* define coefficients going into equation */
/* a(j)*u(j+1)+b(j)*u(j)+c(j)*u(j-1)=d(j) */
    i__1 = *jmax - 1;
    for (j = 2; j <= i__1; ++j) {
	w1 = sqrt(x[j] * x[j + 1]);
	w2 = sqrt(x[j - 1] * x[j]);
/*  w1 is x(j+1/2) */
/*  w2 is x(j-1/2) */
	a[j - 1] = -c20 * c2[j] * (*theta / *deltal / w1 + .5);
	t1 = -c20 * c2[j] * (.5 - *theta / *deltal / w1);
	t2 = c20 * c2[(int) (j - 1.)] * (*theta / *deltal / w2 + .5);
	t3 = pow_dd(&x[j], &c_b11) * (*tautom * bet[j]);
	b[j - 1] = t1 + t2 + t3;
	c__[j - 1] = c20 * c2[j - 1] * (.5 - *theta / *deltal / w2);
	d__[j - 1] = x[j] * dphdot[j];
/* L1: */
    }
/* define constants going into boundary terms */
/* u(1)=aa*u(2) (zero flux at lowest energy) */
/* u(jx2) given from region 2 above */
    x32 = sqrt(x[1] * x[2]);
    aa = (*theta / *deltal / x32 + .5) / (*theta / *deltal / x32 - .5);

/* zero flux at the highest energy */
    u[*jmax - 1] = 0.;

/* invert tridiagonal matrix */
    alp[1] = b[1] + c__[1] * aa;
    gam[1] = a[1] / alp[1];
    i__1 = *jmax - 1;
    for (j = 3; j <= i__1; ++j) {
	alp[j - 1] = b[j - 1] - c__[j - 1] * gam[j - 2];
	gam[j - 1] = a[j - 1] / alp[j - 1];
/* L2: */
    }
    g[1] = d__[1] / alp[1];
    i__1 = *jmax - 2;
    for (j = 3; j <= i__1; ++j) {
	g[j - 1] = (d__[j - 1] - c__[j - 1] * g[j - 2]) / alp[j - 1];
/* L3: */
    }
    g[*jmax - 2] = (d__[*jmax - 2] - a[*jmax - 2] * u[*jmax - 1] - c__[*jmax
	    - 2] * g[*jmax - 3]) / alp[*jmax - 2];
    u[*jmax - 2] = g[*jmax - 2];
    i__1 = *jmax - 1;
    for (j = 3; j <= i__1; ++j) {
	jj = *jmax + 1 - j;
	u[jj - 1] = g[jj - 1] - gam[jj - 1] * u[jj];
/* L4: */
    }
    u[0] = aa * u[1];
/* compute new value of dph(x) and new value of dphesc(x) */
    i__1 = *jmax;
    for (j = 1; j <= i__1; ++j) {
	dphesc[j] = x[j] * x[j] * u[j - 1] * bet[j] * *tautom;
/* L9: */
    }
    return 0;
} /* f_thermlc__ */

static int f_thcompton__(double *tempbb, double *theta,
	double *gamma, double *x, int *jmax, double *sptot)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Builtin functions */
    double pow_dd(double *, double *), sqrt(double), log(
	    double), d_lg10(double *), exp(double);

    /* Local variables */
    static int j;
    static double w;
    static double c2[900], w1, z1, z2, z3, z4, z5, z6, xr, arg, bet[
	    900], rel[900];
    static int jnr;
    static double flz, xnr;
    static int jrel;
    static double xmin, xmax, delta, taukn, deltal, dphesc[900], planck,
	    dphdot[900];
    static int jmaxth;
    static double tautom;

/*  version: January 96 */


/*     Thermal Comptonization; solves Kompaneets eq. with some */
/*     relativistic corrections. See Lightman & Zdziarski (1987), ApJ */
/*     The seed spectrum is blackbody. */
/*  input parameters: */
/* use internally Thomson optical depth */
    /* Parameter adjustments */
    --sptot;
    --x;

    /* Function Body */
    d__1 = *gamma + .5;
    tautom = sqrt(3. / (*theta * (pow_dd(&d__1, &c_b2) - 2.25)) + 2.25) - 1.5;

/* clear arrays (important for repeated calls) */
    for (j = 1; j <= 900; ++j) {
	dphesc[j - 1] = 0.;
	dphdot[j - 1] = 0.;
	rel[j - 1] = 0.;
	bet[j - 1] = 0.;
	c2[j - 1] = 0.;
	sptot[j] = 0.;
/* L31: */
    }

/* JMAX - # OF PHOTON ENERGIES */

/* delta is the 10-log interval of the photon array. */
    delta = .02;
    deltal = delta * log(10.);
    xmin = *tempbb * 1e-4;
    xmax = *theta * 40.;
/* Computing MIN */
    d__1 = xmax / xmin;
    i__1 = 899, i__2 = (int) (d_lg10(&d__1) / delta) + 1;
    *jmax = int_min(i__1,i__2);

/* X - ARRAY FOR PHOTON ENERGIES */

    i__1 = *jmax + 1;
    for (j = 1; j <= i__1; ++j) {
	d__1 = (j - 1) * delta;
	x[j] = xmin * pow_dd(&c_b8, &d__1);
/* L4: */
    }

/* compute c2(x), and rel(x) arrays */
/* c2(x) is the relativistic correction to Kompaneets equation */
/* rel(x) is the Klein-Nishina cross section divided by the Thomson crossection */
    i__1 = *jmax;
    for (j = 1; j <= i__1; ++j) {
	w = x[j];
/* c2 is the Cooper's coefficient calculated at w1 */
/* w1 is x(j+1/2) (x(i) defined up to jmax+1) */
	w1 = sqrt(x[j] * x[j + 1]);
	c2[j - 1] = (double) (pow_dd(&w1, &c_b10) / (w1 * 4.6 + 1. + w1 * 1.1 *
		w1));
	if (w <= .05) {
/* use asymptotic limit for rel(x) for x less than 0.05 */
	    rel[j - 1] = (double) (1. - w * 2. + w * 26. * w / 5.);
	} else {
	    z1 = (w + 1.) / pow_dd(&w, &c_b11);
	    z2 = w * 2. + 1.;
	    z3 = log(z2);
	    z4 = w * 2. * (w + 1.) / z2;
	    z5 = z3 / 2. / w;
	    z6 = (w * 3. + 1.) / z2 / z2;
	    rel[j - 1] = (double) ((z1 * (z4 - z3) + z5 - z6) * .75);
	}
/* L500: */
    }
/* the thermal emision spectrum */
/* Computing MIN */
    d__1 = *tempbb * 50. / xmin;
    i__1 = 900, i__2 = (int) (d_lg10(&d__1) / delta);
    jmaxth = int_min(i__1,i__2);
    if (jmaxth > *jmax) {
/*           print *,'thcomp: ',jmaxth,jmax */
	jmaxth = *jmax;
    }
    d__1 = pi * *tempbb;
    planck = 15. / pow_dd(&d__1, &c_b10);
    i__1 = jmaxth;
    for (j = 1; j <= i__1; ++j) {
	dphdot[j - 1] = planck * pow_dd(&x[j], &c_b2) / (exp(x[j] / *tempbb) 
		- 1.);
/* L5: */
    }

/* compute beta array, the probability of escape per Thomson time. */
/* bet evaluated for spherical geometry and nearly uniform sources. */
/* Between x=0.1 and 1.0, a function flz modifies beta to allow */
/* the increasingly large energy change per scattering to gradually */
/* eliminate spatial diffusion */
    d__1 = .1 / xmin;
    jnr = (int) (d_lg10(&d__1) / delta + 1.);
/* Computing MIN */
    i__1 = jnr, i__2 = *jmax - 1;
    jnr = int_min(i__1,i__2);
    d__1 = 1. / xmin;
    jrel = (int) (d_lg10(&d__1) / delta + 1.);
    jrel = int_min(jrel,*jmax);
    xnr = x[jnr];
    xr = x[jrel];
    i__1 = jnr - 1;
    for (j = 1; j <= i__1; ++j) {
	taukn = tautom * rel[j - 1];
	bet[j - 1] = 1. / tautom / (taukn / 3. + 1.);
/* L501: */
    }
    i__1 = jrel;
    for (j = jnr; j <= i__1; ++j) {
	taukn = tautom * rel[j - 1];
	arg = (x[j] - xnr) / (xr - xnr);
	flz = 1. - arg;
	bet[j - 1] = 1. / tautom / (taukn / 3. * flz + 1.);
/* L600: */
    }
    i__1 = *jmax;
    for (j = jrel + 1; j <= i__1; ++j) {
	bet[j - 1] = 1. / tautom;
/* L601: */
    }

    f_thermlc__(&tautom, theta, &deltal, &x[1], jmax, dphesc, dphdot, bet, c2)
	    ;

/*     the spectrum in E F_E */
    i__1 = *jmax - 1;
    for (j = 1; j <= i__1; ++j) {
	sptot[j] = dphesc[j - 1] * pow_dd(&x[j], &c_b2);
/*          write(1,*) x(j), sptot(j) */
/* L497: */
    }
/*     the input spectrum */
/*      do 498 j=1,jmaxth */
/*         write(2,*) x(j), dphdot(j)*x(j)**2 */
/* 498  continue */

    return 0;
} /* f_thcompton__ */


static int f_thdscompton__(double *tempbb, double *theta,
	double *gamma, double *x, int *jmax, double *sptot)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Builtin functions */
    double pow_dd(double *, double *), sqrt(double), log(
	    double), d_lg10(double *);

    /* Local variables */
    static int j;
    static double w;
    extern /* Subroutine */ int f_thermlc__(double *, double *,
	    double *, double *, int *, double *, double *,
	     double *, double *);
    static double c2[900], w1, z1, z2, z3, z4, z5, z6;
    static int ne;
    static double xr, ear[5001], arg, bet[900];
    static int ifl;
    static double rel[900];
    static int jnr;
    static double flz, xnr;
    static int jrel;
    static double xmin, xmax, delta, parth[10], taukn, deltal, dphesc[900]
	    , dphdot[900];
    static int jmaxth;
    static double photar[5000], photer[5000], tautom;

/*  version: January 96 */


/*     Thermal Comptonization; solves Kompaneets eq. with some */
/*     relativistic corrections. See Lightman & Zdziarski (1987), ApJ */
/*     The seed spectrum is DISK blackbody. */
/*  input parameters: */
/* use internally Thomson optical depth */
    /* Parameter adjustments */
    --sptot;
    --x;

    /* Function Body */
    d__1 = *gamma + .5;
    tautom = sqrt(3. / (*theta * (pow_dd(&d__1, &c_b2) - 2.25)) + 2.25) - 1.5;

/* clear arrays (important for repeated calls) */
    for (j = 1; j <= 900; ++j) {
	dphesc[j - 1] = 0.;
	dphdot[j - 1] = 0.;
	rel[j - 1] = 0.;
	bet[j - 1] = 0.;
	c2[j - 1] = 0.;
	sptot[j] = 0.;
/* L31: */
    }

/* JMAX - # OF PHOTON ENERGIES */

/* delta is the 10-log interval of the photon array. */
    delta = .02;
    deltal = delta * log(10.);
    xmin = *tempbb * 1e-4;
    xmax = *theta * 40.;
/* Computing MIN */
    d__1 = xmax / xmin;
    i__1 = 899, i__2 = (int) (d_lg10(&d__1) / delta) + 1;
    *jmax = int_min(i__1,i__2);

/* X - ARRAY FOR PHOTON ENERGIES */

    i__1 = *jmax + 1;
    for (j = 1; j <= i__1; ++j) {
	d__1 = (j - 1) * delta;
	x[j] = xmin * pow_dd(&c_b8, &d__1);
/* L4: */
    }

/* compute c2(x), and rel(x) arrays */
/* c2(x) is the relativistic correction to Kompaneets equation */
/* rel(x) is the Klein-Nishina cross section divided by the Thomson crossection */
    i__1 = *jmax;
    for (j = 1; j <= i__1; ++j) {
	w = x[j];
/* c2 is the Cooper's coefficient calculated at w1 */
/* w1 is x(j+1/2) (x(i) defined up to jmax+1) */
	w1 = sqrt(x[j] * x[j + 1]);
	c2[j - 1] = (double) (pow_dd(&w1, &c_b10) / (w1 * 4.6 + 1. + w1 * 1.1 *
		w1));
	if (w <= .05) {
/* use asymptotic limit for rel(x) for x less than 0.05 */
	    rel[j - 1] = (double) (1 - w * 2 + w * 26 * w / 5);
	} else {
/* Computing 3rd power */
	    d__1 = w;
	    z1 = (w + 1) / (d__1 * (d__1 * d__1));
	    z2 = w * 2 + 1;
	    z3 = log(z2);
	    z4 = w * 2 * (w + 1) / z2;
	    z5 = z3 / 2 / w;
	    z6 = (w * 3 + 1) / z2 / z2;
	    rel[j - 1] = (double) ((z1 * (z4 - z3) + z5 - z6) * .75);
	}
/* L500: */
    }
/* the thermal emision spectrum */
/* Computing MIN */
    d__1 = *tempbb * 50. / xmin;
    i__1 = 900, i__2 = (int) (d_lg10(&d__1) / delta);
    jmaxth = int_min(i__1,i__2);
    if (jmaxth > *jmax) {
/*           print *,'thcomp: ',jmaxth,jmax */
	jmaxth = *jmax;
    }
/*        planck=15/(pi*tempbb)**4 */
/*        do 5 j=1,jmaxth */
/*          dphdot(j)=planck*x(j)**2/(exp(x(j)/tempbb)-1) */
/*  5     continue */
    i__1 = jmaxth - 1;
    for (j = 1; j <= i__1; ++j) {
	ear[j - 1] = sqrt(x[j] * x[j + 1]) * 511.;
    }
    parth[0] = *tempbb * 511.;
    ne = jmaxth - 2;
    f_xsdskb__(ear, &ne, parth, &ifl, photar, photer);
    i__1 = ne;
    for (j = 1; j <= i__1; ++j) {
	dphdot[j] = photar[j - 1] * 511. / (ear[j] - ear[j - 1]);
    }
    jmaxth = ne + 1;
    dphdot[0] = dphdot[1];

/* compute beta array, the probability of escape per Thomson time. */
/* bet evaluated for spherical geometry and nearly uniform sources. */
/* Between x=0.1 and 1.0, a function flz modifies beta to allow */
/* the increasingly large energy change per scattering to gradually */
/* eliminate spatial diffusion */
    d__1 = .1 / xmin;
    jnr = (int) (d_lg10(&d__1) / delta + 1);
/* Computing MIN */
    i__1 = jnr, i__2 = *jmax - 1;
    jnr = int_min(i__1,i__2);
    d__1 = 1. / xmin;
    jrel = (int) (d_lg10(&d__1) / delta + 1);
    jrel = int_min(jrel,*jmax);
    xnr = x[jnr];
    xr = x[jrel];
    i__1 = jnr - 1;
    for (j = 1; j <= i__1; ++j) {
	taukn = tautom * rel[j - 1];
	bet[j - 1] = 1 / tautom / (taukn / 3 + 1);
/* L501: */
    }
    i__1 = jrel;
    for (j = jnr; j <= i__1; ++j) {
	taukn = tautom * rel[j - 1];
	arg = (x[j] - xnr) / (xr - xnr);
	flz = 1 - arg;
	bet[j - 1] = 1 / tautom / (taukn / 3 * flz + 1);
/* L600: */
    }
    i__1 = *jmax;
    for (j = jrel + 1; j <= i__1; ++j) {
	bet[j - 1] = 1 / tautom;
/* L601: */
    }

    f_thermlc__(&tautom, theta, &deltal, &x[1], jmax, dphesc, dphdot, bet, c2)
	    ;

/*     the spectrum in E F_E */
    i__1 = *jmax - 1;
    for (j = 1; j <= i__1; ++j) {
/* Computing 2nd power */
	d__1 = x[j];
	sptot[j] = dphesc[j - 1] * (d__1 * d__1);
/*          write(1,*) x(j), sptot(j) */
/* L497: */
    }
/*      print *,'jmax: ',jmax,jmaxth */
/*      open(33,file='spec.dat') */
/* c     the input spectrum */
/*      do 498 j=1,min(jmaxth,jmax-1) */
/*         write(33,*) 511*x(j), dphdot(j)*x(j), dphesc(j)*x(j) */
/* 498  continue */
/*      close(33) */
    return 0;
} /* f_thdscompton__ */

/* ------------------------------------------------------------------ c */
static double f_spp__(double *y, double *xnonth, int *nnonth,
	double *spnth)
{
    /* Initialized data */

    static int ih = 2;

    /* System generated locals */
    double ret_val;

    /* Local variables */
    static int il;
    static double xx;

    /* Parameter adjustments */
    --spnth;
    --xnonth;

    /* Function Body */
    xx = 1 / *y;
    if (xx < xnonth[ih]) {
	ih = 2;
    }
    while(ih < *nnonth && xx > xnonth[ih]) {
	++ih;
    }
    il = ih - 1;
    ret_val = spnth[il] + (spnth[ih] - spnth[il]) * (xx - xnonth[il]) / (
	    xnonth[ih] - xnonth[il]);
    return ret_val;
} /* f_spp__ */

void c_donthcomp(double *ear, int ne, double* param, double *photar){
    /* Initialized data */
	double prim[ne];

    static double pa0[5] = { 9999.,9999.,9999.,9999.,9999. };

    /* System generated locals */
    int i__1;
    double d__1, d__2;

    /* Builtin functions */
    double pow_dd(double *, double *);

    /* Local variables */
    static int i__, j, n, jl, np;
    static double xn;
    static int fl0;
    static int nth;
    static double xth[900], spt[900];
    static double z_red__;
    static double normfac, normlum;


/*     driver for the Comptonization code solving Kompaneets equation */
/*     seed photons - (disc) blackbody */
/*     reflection + Fe line with smearing */

/*     Version optimized for a number of data files but with the same values */
/*     of parameters: */
/*        computes the model when ifl=1, in broad energy range, and rebins */
/*        it only for subsequent files */

/*     number of model parameters: 16 */
/*     1: photon spectral index */
/*     2: plasma temperature in keV */
/*     3: (disc)blackbody temperature in keV */
/*     4: type of seed spectrum (0-blackbody, 1-diskbb) */
/*     5: redshift */
/*      INCLUDE 'xspec.inc' */
/*     SPP is external (xansrc/functions/xsnteea.f) */
    /* Parameter adjustments */
    --photar;
    --param;

    /* Function Body */

/*     xtot is the energy array (units m_e c^2) */
/*     spnth is the nonthermal spectrum alone (E F_E) */
/*     sptot is the total spectrum array (E F_E), = spref if no reflection */
    z_red__ = param[5];
/*  calculate internal source spectrum if necessary */
/*  (only for the first file) */
    np = 5;
    i__1 = np;
    for (n = 1; n <= i__1; ++n) {
	if (param[n] != pa0[n - 1]) {
	    fl0 = 1;
	}
    }
    if (fl0) {
	if (param[4] < .5) {
	    d__1 = param[3] / 511.;
	    d__2 = param[2] / 511.;
	    f_thcompton__(&d__1, &d__2, &param[1], xth, &nth, spt);
	} else {
	    d__1 = param[3] / 511.;
	    d__2 = param[2] / 511.;
	    f_thdscompton__(&d__1, &d__2, &param[1], xth, &nth, spt);
	}
    }
    xn = (z_red__ + 1) / 511.;
    d__1 = 1 / xn;
    normfac = 1 / f_spp__(&d__1, xth, &nth, spt);
/* Calculate luminosity normalization  (used in another model!) */
    normlum = 0.;
    i__1 = nth - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	normlum += (spt[i__ - 1] / xth[i__ - 1] + spt[i__ - 2] / xth[i__ - 2])
		 * .5 * (xth[i__ - 1] - xth[i__ - 2]);
    }
    normlum *= normfac;
/*     zero arrays */
    i__1 = ne;
    for (i__ = 1; i__ <= i__1; ++i__) {
	photar[i__] = 0.;
	prim[i__] = 0.;
    }
    prim[0] = 0.;

/*     put primary into final array only if scale >= 0. */
    j = 1;
    i__1 = ne;
    for (i__ = 0; i__ <= i__1; ++i__) {
	while(j <= nth && xth[j - 1] * 511. < ear[i__] * (z_red__ + 1)) {
	    ++j;
	}
	if (j <= nth) {
	    if (j > 1) {
		jl = j - 1;
		prim[i__] = spt[jl - 1] + (ear[i__] / 511. * (z_red__ + 1) -
			xth[jl - 1]) * (spt[jl] - spt[jl - 1]) / (xth[jl] -
			xth[jl - 1]);
	    } else {
		prim[i__] = spt[0];
	    }
	}
    }
    i__1 = ne;
    for (i__ = 1; i__ <= i__1; ++i__) {
	photar[i__] = (prim[i__] / pow_dd(&ear[i__], &c_b2) + prim[i__ - 1] /
		pow_dd(&ear[i__ - 1], &c_b2)) * .5 * (ear[i__] - ear[i__ - 1])
		 * normfac;
    }
    return;
} /* f_donthcomp__ */
