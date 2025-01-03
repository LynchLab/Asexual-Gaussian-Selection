/* 

This program estimates the features of a multiple-allele locus in a finite population, allowing for exponential or quadratic (or any other exponent) selection 
and reversible mutation. The goal is to estimate the steady-state distribution of alternative allele types over time. The equilibrium long-term result is a balance between 
the forces of mutation, selection, and random genetic drift. 

***** This version has three types of sites simultaneously present denoted 1, 2, and 3.

***** This version allows the use of arbitrary numbers of sites and stength of selection, and also an arbitrary fitness function of the form exp(-s * (n^powerexp)).
      Powerexp = 1.0 is the negative exponential model; Powerexp = 2.0 is the half-Gaussian model; but any other exponent is allowable. 

The actual population size is assumed to be equal to the effective size (N, the number of new individuals drawn each generation), 
which is constant in time, although the genetic effective size influenced by linkage naturally evolves.

The population experiences sequential episodes of mutation, selection, and random genetic drit.

Haploidy is assumed, and there is no recombination.

Allele designation:	the general code applies to any trait with a series of sites, each with biallelic states -/+.
	-----....-----, for the most deleterious allele.
	+++++....+++++, for the most fit allele.

But because there is complete linkage, order does not matter, and a fraction of sites is allocated to the minor vs. major component.

The population is regularly censused after reproduction, at which point the population has size N.

All mutation and selection processes are treated deterministically, prior to random sampling.

Mutation:
	There are just two types of mutations at each type of site: - to +, u01 (beneficial); and + to -, u10 (deleterious); the same rates are assumed at all loci.
	
Selection:
	Fitness is determined by a function that has to be set internally. 

The run starts with an allele-frequency distribution that can be set near the top of the program.

After a burnin, statistics are then taken at intervals of Ne/xinc generations, with a total of ngen sampling intervals.

The rates of mutation and strength of selection are defined by input parameters entered immediately below.

NOTE THAT NEAR THE TOP OF THE PROGRAM, A SCALING FACTOR (kfac) CAN BE UTILIZED TO SCALE UP THE SELECTION COEFFICIENTS. IF THE NE IS SCALED DOWN AND THE 
MUTATION RATES UP BY THE SAME FACTOR, THIS ENABLES THE OVERALL PROGRAM TO RUN FASTER, AS ALL IS A FUNCTION OF THE PRODUCTS NE*S AND NE*U.

********** CARE NEEDS TO BE TAKEN IN SETTING kfac FOR EACH POPULATION SIZE TO BE SURE THAT: 
	1) THE MUTATION RATE * NUMBER OF SITES * kfac IS SMALLER THAN 1.0, IDEALLY < 0.1 (TO AVOID MUTATION TRANSITION PROBABILITIES EXCEEDING 1.0); AND 
	2) THE SELECTION COEFFICIENT * kfac IS SMALLER THAN 1.0 (TO AVOID NEGATIVE FITNESSES); THIS IS TOUCHY BECAUSE THE ACTUAL STRENGTH OF SELECTION 
	   DEPENDS ON THE LOCATION OF THE PHENOTYPE WHEN THE MODEL IS ANYTHING OTHER THAN THE EXPONENTIAL.

********** THERE ARE THREE PLACES WHERE THE NAME OF THE OUTPUT FILE DATAOUT... NEEDS TO BE EDITTED FOR THE PRINTOUTS.
	

A SERIES OF 21 COMBINATIONS OF POPULATION SIZE (efpopn) AND MUTATION RATE (delmutrate) (USING SCALING FROM LYNCH ET AL. 2016), COVERING THE FULL RANGE OF 
NATURAL POPULATION VALUES, IS GIVEN INTERNALLY. THE LOOP IMMEDIATELY BELOW THIS (itera) NEEDS TO BE SET TO DETERMINE THE RANGE OF POPULATION SIZES TO RUN. 
 



 ********************************************************************************************************************** */

#define mutbeta		1.0						/* ratio of beneficial to deleterious mutation rates */

#define ell1		3333						/* total number of sites for the trait with small effects */
#define ell2		333						/* total number of sites for the trait with medium effects */
#define ell3		33
#define sco1		0.00001						/* selection coefficient for small-effect alleles */
#define sco2		0.0001 				    /* selection coefficient for medium-effect alleles */
#define sco3		0.001				/* selection coefficient for major-effect alleles. BE SURE TO PLACE IN INCREASING ORDER */

#define powerexp	2.0						/* 1.0 means an exponential function; 2.0 means a half-Gaussian. */

/* #define thetaopt		((double) ell3) + (((double) ell2) * (sco2/sco3)) + (((double) ell1) * (sco1/sco3)) */		/* THIS SETTING APPLIES TO THE HALF-GAUSSIAN */

#define thetaopt    50.0                                                                       /*  THIS SETTING APPLIES TO THE INTERMEDIATE OPTIMUM */

#define xinc		10						/* statistics to be recorded every ne/xinc generations */

#define burnin		1000					/* number of initial burn-in sampling increments, ignored in the statistics */

#define tintprint	1000					/* number of sampling increments between screen printing */



#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include	<time.h>
#include    <gsl/gsl_rng.h>
#include    <gsl/gsl_randist.h>
#include    <string.h>


/* ********************************************************************************************************** */



/* point to the output file */

FILE *stream;
#define MAX_FILENAME_LENGTH 100
char filename[MAX_FILENAME_LENGTH];

int main(int argc, char *argv[])
{

int f0, f1;
if (argc > 1)
{
    f0 = atoi(argv[1]);
    f1=f0; 
    snprintf(filename, MAX_FILENAME_LENGTH, "dataout2006_%d.txt", f0); 
}
else
{
    f0 = 1;
    f1 = 21;
    snprintf(filename, MAX_FILENAME_LENGTH, "dataout2006.txt");
}

static gsl_rng* rand_new;                                                       
static gsl_rng_type * T;                                                        
gsl_rng_env_setup();                                                            
if(!rand_new)                                                                   
{                                                                               
    rand_new = gsl_rng_alloc(gsl_rng_taus2);                                    
    gsl_rng_set(rand_new, time(NULL));                                          
}     

                                                                                
 
   

/* ***************************************************************************************** */



/* MAIN BODY OF PROGRAM. */


int ig, ig1, ig2, ig3, jg;								/* counters for the classes, with 0 meaning no + alleles, and ell meaning all + alleles */

double start, stop, time;

long igen;												/* generation counter */
double tint;											/* interval number for statistics */
long increment;											/* increment between surveys */
long tcount;											/* counter for printing to screen to monitor simulation progress */
long counter;											/* initial number of surveys to be skipped */

double ngens;                                           /* time iterations in run */

double phi;									        	/* ratio of beneficial and deleterious mutation rates per site */

int ellfac;                                             /* highest number of sites for the three types, for obtaining the normalizing factor */
double scofac;                                          /* highest selection coefficient for the three types, for obtaining the normalizing factor */
double maxselco;										/* maximum strength of selection across adjacent classes, used to get kfac */

double selco1, selco2, selco3;							/* selection coefficients associated with major- and minor-effect loci */

double exterm, pseq1, pseq2, pseq3;					    /* expected frequencies under the sequential model (assumes no interference) */

long efpopn[40];										/* effective population sizes to run -- NOTE THESE ARE SCALED DOWN ACCORDING TO kfac TO INCREASE RUN SPEED */
double delmutrate[40];								    /* deleterious mutation rates to run -- NOTE THESE ARE SCALED UP ACCORDING TO kfac TO INCREASE RUN SPEED */
int scalef[40];											/* SCALING FACTORS FOR SPEEDING UP RUNS */
double rlng[40];										/* number of iterations to run the program */

long ne;												/* effective population size, from efpopn[] */
double u10, u01;									    /* deleterious and beneficial mutation rates, from delmutrate[] */
int kfac;												/* scaling factor for speeding up runs, from scalef[] */

double kfacn, kfacs, kfacu;                             /* possible scaling factors based on population size, selection coefficients, and mutation rates */
    
int itera;												/* counter for population-size/mutation-rate iterations */

double dif;												/* measure of phenotype scaled to modified thetaopt */

double meanfit;									        /* mean fitness */

double mut01, mutdown01, mute1, mutupe1;				/* parameters for mutational transition frequencies between classes */
double mut02, mutdown02, mute2, mutupe2;
double mut03, mutdown03, mute3, mutupe3;

double totw, totwsq;								    /* summations for grand means and variances of fitness */

double totmin0, meanminhet, neuthet;				    /* statistics for computing the long-term mean heterozygosity at a single neutral minor site -- ***** ell1 must be set to 1, and selco1 = 0.0 to be meaningful ***** */

double totig1, totigsq1, totwvarig1;					/* summations for grand means, variances, and higher-order within-population moments */
double totig2, totigsq2, totwvarig2;
double totig3, totigsq3, totwvarig3;
double totw3mom1, totw4mom1, totw3mom2, totw4mom2, totw3mom3, totw4mom3;

int newlow1, newhigh1, newlow2, newhigh2, newlow3, newhigh3;		/* running settings for upper and lower genoytpic states for allele counts for the three site types */
int lowig1, highig1, lowig2, highig2, lowig3, highig3;
int lowmut1, highmut1, lowmut2, highmut2, lowmut3, highmut3;

double sump;										    /* sum of frequencies */

double pp;											    /* probability associated with the binomial for drift */
long ntot;												/* integer associated with the binomial for drift */
long draw;												/* drift drawn from the binomial */

double meanw, meanig1, meanig2, meanig3;		        /* generational means for fitness and allelic class */

double ssqig1, wvarig1;									/* sum of squares and within-generation variance of trait */
double ssqig2, wvarig2;
double ssqig3, wvarig3;

double grandmeanig1, varig1;							/* mean and variances of numbers of minor and major alleles*/
double grandmeanig2, varig2;
double grandmeanig3, varig3;

double grandmeanw, varw;							    /* mean and variance of fitness */

double meanperform;								        /* mean relative performance under additive model */

double neestimate;

double fracpine;                                        /* effective population size based on pis relative to input ne */ 

double mean1, mean2, mean3;                             /* grand average frequencies for minor and major alleles */

double bne1, bne2, bne3;                                /* Ne back-estimated from the Bulmer sequential model equation, using observed mean + allele frequencies */
double fracbne1, fracbne2, fracbne3;

double bmean, bmnsq, bvar, bsd1, bsd2, bsd3;            /* statistics for the between-generation variance of mean numbers of + alleles */

double perform, totperform, totperformsq, bsdperform;					/* statistics for the mean and sd among-generations in performance */
double devi, totdevi, totdevisq, meandevi, bsddevi;						/* statistics for the mean and sd among-generations of deviation of performance from optimum */
double absdevi, totabsdevi, totabsdevisq, meanabsdevi, bsdabsdevi;		/* statistics for the mean and sd among-generations of absolute deviation of performance from optimum */

double w3momig1, w3momig2, w3momig3;					/* third-order within-generation moment */
double w4momig1, w4momig2, w4momig3;					/* fourth-order within-generation moment */



/* Format the genotypic arrys. */

double *mutransup1 = (double*)malloc((ell1+2) * sizeof(double));
double *mutransdown1 = (double*)malloc((ell1+2) * sizeof(double));
double *mutrans01 = (double*)malloc((ell1+2) * sizeof(double));

double *mutransup2 = (double*)malloc((ell2 + 2) * sizeof(double));
double *mutransdown2 = (double*)malloc((ell2 + 2) * sizeof(double));
double *mutrans02 = (double*)malloc((ell2 + 2) * sizeof(double));

double *mutransup3 = (double*)malloc((ell3 + 2) * sizeof(double));
double *mutransdown3 = (double*)malloc((ell3 + 2) * sizeof(double));
double *mutrans03 = (double*)malloc((ell3 + 2) * sizeof(double));

double *sumfreq1 = (double*)malloc((ell1 + 2) * sizeof(double)); 
double *sumfreq2 = (double*)malloc((ell2 + 2) * sizeof(double)); 
double *sumfreq3 = (double*)malloc((ell3 + 2) * sizeof(double));

double *meanpfreq1 = (double*)malloc((ell1 + 2) * sizeof(double));
double *meanpfreq2 = (double*)malloc((ell2 + 2) * sizeof(double)); 
double *meanpfreq3 = (double*)malloc((ell3 + 2) * sizeof(double));



double ***relw = (double ***)malloc((ell1 + 2) * (sizeof(double **)));
	for (ig = 0; ig < (ell1 + 2); ig++) {
		relw[ig] = (double **)malloc((ell2 + 2) * (sizeof (double *)));
	    for (jg = 0; jg < (ell2 + 2); jg++) {
			relw[ig][jg] = (double *)malloc((ell3 + 2) * (sizeof (double))); }	}

double ***corwfit = (double ***)malloc((ell1 + 2) * (sizeof(double **)));
	for (ig = 0; ig < (ell1 + 2); ig++) {
		corwfit[ig] = (double **)malloc((ell2 + 2) * (sizeof (double *)));
	    for (jg = 0; jg < (ell2 + 2); jg++) {
		corwfit[ig][jg] = (double *)malloc((ell3 + 2) * (sizeof (double))); }	}

double ***p0 = (double ***)malloc((ell1 + 2) * (sizeof(double **)));
	for (ig = 0; ig < (ell1 + 2); ig++) {
		p0[ig] = (double **)malloc((ell2 + 2) * (sizeof (double *)));
	    for (jg = 0; jg < (ell2 + 2); jg++) {
			p0[ig][jg] = (double *)malloc((ell3 + 2) * (sizeof (double))); 	}}

double ***pmutm1 = (double ***)malloc((ell1 + 2) * (sizeof(double **)));
	for (ig = 0; ig < (ell1 + 2); ig++) {
		pmutm1[ig] = (double **)malloc((ell2 + 2) * (sizeof (double *)));
	    for (jg = 0; jg < (ell2 + 2); jg++) {
			pmutm1[ig][jg] = (double *)malloc((ell3 + 2) * (sizeof (double)));} }

double ***pmutm2 = (double ***)malloc((ell1 + 2) * (sizeof(double **)));
	for (ig = 0; ig < (ell1 + 2); ig++) {
		pmutm2[ig] = (double **)malloc((ell2 + 2) * (sizeof (double *)));
	    for (jg = 0; jg < (ell2 + 2); jg++) {
			pmutm2[ig][jg] = (double *)malloc((ell3 + 2) * (sizeof (double))); }}	

double ***pmutm = (double ***)malloc((ell1 + 2) * (sizeof(double **)));
	for (ig = 0; ig < (ell1 + 2); ig++) {
		pmutm[ig] = (double **)malloc((ell2 + 2) * (sizeof (double *)));
	    for (jg = 0; jg < (ell2 + 2); jg++) {
			pmutm[ig][jg] = (double *)malloc((ell3 + 2) * (sizeof (double))); }}

double ***psel = (double ***)malloc((ell1 + 2) * (sizeof(double **)));
	for (ig = 0; ig < (ell1 + 2); ig++) {
		psel[ig] = (double **)malloc((ell2 + 2) * (sizeof (double *)));
	    for (jg = 0; jg < (ell2 + 2); jg++) {
			psel[ig][jg] = (double *)malloc((ell3 + 2) * (sizeof (double))); 	}}

double ***pgtypexp = (double ***)malloc((ell1 + 2) * (sizeof(double **)));
	for (ig = 0; ig < (ell1 + 2); ig++) {
		pgtypexp[ig] = (double **)malloc((ell2 + 2) * (sizeof (double *)));
	    for (jg = 0; jg < (ell2 + 2); jg++) {
			pgtypexp[ig][jg] = (double *)malloc((ell3 + 2) * (sizeof (double))); }	}




/* Open the output file. */

remove("dataout2006.txt ");


/* Effective population sizes to run */

efpopn[21] = 205636404;
efpopn[20] = 157180983;
efpopn[19] = 117462706;
efpopn[18] = 85822271;
efpopn[17] = 61305579;
efpopn[16] = 42815399;
efpopn[15] = 29234791;
efpopn[14] = 19516413;
efpopn[13] = 12737963;
efpopn[12] = 8128305;
efpopn[11] = 5071075;
efpopn[10] = 3093000;
efpopn[9] = 1844591;
efpopn[8] = 1075474;
efpopn[7] = 613056;
efpopn[6] = 341665;
efpopn[5] = 186166;
efpopn[4] = 99174;
efpopn[3] = 51654;
efpopn[2] = 26300;
efpopn[1] = 13095;


/* Associated mutation rates */

delmutrate[21] = 0.0000000005290;
delmutrate[20] = 0.0000000006488;
delmutrate[19] = 0.0000000008096;
delmutrate[18] = 0.000000001028;
delmutrate[17] = 0.000000001327;
delmutrate[16] = 0.000000001743;
delmutrate[15] = 0.000000002330;
delmutrate[14] = 0.000000003167;
delmutrate[13] = 0.000000004380;
delmutrate[12] = 0.000000006163;
delmutrate[11] = 0.000000008821;
delmutrate[10] = 0.00000001284;
delmutrate[9] = 0.00000001900;
delmutrate[8] = 0.00000002870;
delmutrate[7] = 0.00000004390;
delmutrate[6] = 0.00000006850;
delmutrate[5] = 0.0000001090;
delmutrate[4] = 0.0000001750;
delmutrate[3] = 0.0000002880;
delmutrate[2] = 0.0000004810;
delmutrate[1] = 0.0000008170; 


/* Number of sampling increments in run; each increment is (ne/10) generations */

rlng[21] = 500.0;
rlng[20] = 500.0;
rlng[19] = 500.0;
rlng[18] = 2000.0;
rlng[17] = 2000.0;
rlng[16] = 2000.0;
rlng[15] = 2000.0;
rlng[14] = 2000.0;
rlng[13] = 10000.0;
rlng[12] = 10000.0;
rlng[11] = 10000.0;
rlng[10] = 10000.0;
rlng[9] = 20000.0;
rlng[8] = 40000.0;
rlng[7] = 40000.0;
rlng[6] = 60000.0;
rlng[5] = 600000.0;
rlng[4] = 600000.0;
rlng[3] = 1000000.0;
rlng[2] = 1000000.0;
rlng[1] = 1000000.0;




for (itera = f0; itera <= f1; ++itera) {							/* Start iterations over the set of population sizes and mutation rates. */

stream=fopen(filename, "a");		



/* Set the run length, and determine scaling factor kfac for N, u, and s to speed up run. */
        
    ngens = rlng[itera];

	ne = efpopn[itera];											/* effective population size */

	u10 = delmutrate[itera];									/* deleterious and beneficial mutation rates */
	u01 = u10 * mutbeta;

    kfacn = ((double) ne) / 1000.0;                             /* determine scaling factor kfac for N, such that the actual N run is always > 1000 (below) */
    

	scofac = 0.0;												/* get the potential scaling factor based on selection coefficients, using the largest s, 
																     so that scaled s < 0.1 (below) */

	for (ig3 = 1; ig3 <= ell3; ++ig3) {							/* differences in fitness between adjacent states at site type 3, the major locus type */
																/* THIS NEEDS TO BE EDITTED FOR OTHER FITNESS FUNCTIONS */

		maxselco = exp(-sco3 * pow((thetaopt + 1.0 - (double) ig3), powerexp)) - exp(-sco3 * pow((thetaopt - (double) ig3), powerexp));
		maxselco = abs(maxselco);

		if (maxselco > scofac) {
			scofac = maxselco; 	}
	}
    
    kfacs = 0.1 / scofac;
    
    
    if ((ell2 >= ell1) && (ell2 >= ell3)) {    /* get the potential scaling factor based on the mutation rate, using the largest number of sites, so the max mutation rate < 0.1 */
        ellfac = ell2;  }
    else if ((ell3 >= ell1) && (ell3 >= ell2)) {
        ellfac = ell3;  }
    else {ellfac = ell1;}
    
    kfacu = 0.1 / (u10 * ellfac);

    if (kfacn < 1.0) {
        kfacn = 1.0; }
    if (kfacs < 1.0) {
        kfacs = 1.0; }
    if (kfacu < 1.0) {
        kfacu = 1.0; }

    if ((kfacn < kfacs) && (kfacn < kfacu)) {                   /* of the three potential scaling factors, choose the smallest one */
        kfac = int(kfacn); }
    if ((kfacs < kfacn) && (kfacs < kfacu)) {
        kfac = int(kfacs); }
    if ((kfacu < kfacs) && (kfacu < kfacn)) {
        kfac = int(kfacu); }


	ne = ne / kfac;                                             /* rescale the population sizes and mutation rates */
	u10 = ((double)kfac) * u10;
	u01 = ((double)kfac) * u01;


	selco1 = ((double)kfac) * sco1;                             /* rescale the selection coefficients */
	selco2 = ((double)kfac) * sco2;
	selco3 = ((double)kfac) * sco3;



/* Sequential model expectations for the equilibrium allele frequencies (with no interference), based on the Li-Bulmer model */
	/* NOTE THAT THESE ESTIMATES ONLY MAKE SENSE FOR THE EXPONENTIAL MODEL, FOR OTHERWISE S CHANGES WITH DISTANCE FROM THE OPTIMUM */

	exterm = 2.0 * ((double)ne) * selco1;						
	exterm = mutbeta * exp(exterm);								
	pseq1 = exterm / (1.0 + exterm);
	
	exterm = 2.0 * ((double)ne) * selco2;
	exterm = mutbeta * exp(exterm);
	pseq2 = exterm / (1.0 + exterm);

	exterm = 2.0 * ((double)ne) * selco3;
	exterm = mutbeta * exp(exterm);
	pseq3 = exterm / (1.0 + exterm);



	/* Set the genotypic fitnesses */

	for (ig1 = 0; ig1 <= ell1; ++ig1) {					
		sumfreq1[ig1] = 0.0;                            /* zero the frequency counters for later use */

		for (ig2 = 0; ig2 <= ell2; ++ig2) {
			sumfreq2[ig2] = 0.0;

			for (ig3 = 0; ig3 <= ell3; ++ig3) {
				sumfreq3[ig3] = 0.0;

				dif = ((double)ig3) + (((double)ig2) * (selco2 / selco3)) + (((double)ig1) * (selco1 / selco3));
				
				relw[ig1][ig2][ig3] = exp(-selco3 * pow(thetaopt - dif, powerexp)) ;
							
			    corwfit[ig1][ig2][ig3] = exp(-(selco3 / ((double)kfac))  * pow(thetaopt - dif, powerexp)) ;

			}
		}
	}                                             /* above are fitnesses on the new scale (wfit and relw) and original scale (corw) */ 




	/* Set the mutation constants. */

	phi = u01 / u10;										/* ratio of beneficial and deleterious mutation rates per site */
	 
	mut01 = 1.0 - (((double)ell1) * u01);				    /* fraction remaining in 1 class 0 */
	mutdown01 = u10;										/* fraction of 1 class 1 degrading to 1 class 0 */
	mute1 = 1.0 - (((double)ell1) * u10);				    /* fraction remaining in 1 class ell1 */
	mutupe1 = u01;											/* fraction of 1 class (ell1-1) moving to 1 class ell1 */

	mut02 = 1.0 - (((double)ell2) * u01);				    /* fraction remaining in 2 class 0 */
	mutdown02 = u10;										/* fraction of 2 class 1 degrading to 2 class 0 */
	mute2 = 1.0 - (((double)ell2) * u10);				    /* fraction remaining in 2 class ell2 */
	mutupe2 = u01;											/* fraction of 2 class (ell2-1) moving to 2 class ell2 */

	mut03 = 1.0 - (((double)ell3) * u01);				    /* fraction remaining in 3 class 0 */
	mutdown03 = u10;										/* fraction of 3 class 1 degrading to 3 class 0 */
	mute3 = 1.0 - (((double)ell3) * u10);				    /* fraction remaining in 3 class ell3 */
	mutupe3 = u01;											/* fraction of 3 class (ell3-1) moving to 3 class ell3 */

	for (ig1 = 1; ig1 <= (ell1 - 1); ++ig1) {
		mutransup1[ig1] = u01 * (((double)ell1) - ((double)(ig1 - 1))); 											/* gain of a + from any of the - in next lowest class */
		mutransdown1[ig1] = u10 * ((double)(ig1 + 1));																/* loss of a + from next highest class */
		mutrans01[ig1] = 1.0 - (u10 * ((double)ig1)) - (u01 * (((double)ell1) - ((double)ig1))); 	}		        /* stays unchanged */

	for (ig2 = 1; ig2 <= (ell2 - 1); ++ig2) {
		mutransup2[ig2] = u01 * (((double)ell2) - ((double)(ig2 - 1))); 											
		mutransdown2[ig2] = u10 * ((double)(ig2 + 1));																
		mutrans02[ig2] = 1.0 - (u10 * ((double)ig2)) - (u01 * (((double)ell2) - ((double)ig2))); 	}		        

	for (ig3 = 1; ig3 <= (ell3 - 1); ++ig3) {
		mutransup3[ig3] = u01 * (((double)ell3) - ((double)(ig3 - 1)));
		mutransdown3[ig3] = u10 * ((double)(ig3 + 1));
		mutrans03[ig3] = 1.0 - (u10 * ((double)ig3)) - (u01 * (((double)ell3) - ((double)ig3)));	}




	/* Set the initial genotype frequencies. */

	for (ig1 = 0; ig1 <= ell1; ++ig1) {
		for (ig2 = 0; ig2 <= ell2; ++ig2) {
			for (ig3 = 0; ig3 <= ell3; ++ig3) {
				p0[ig1][ig2][ig3] = 0.0; }}}

	p0[1500][150][15] = 1.0;					/* INITIALIZE STARTING (PRE-BURNIN) FREQUENCIES + */




	/* Initiate the allele frequencies, counters, and test statistics. */

	for (ig1 = 0; ig1 <= ell1; ++ig1) {
		for (ig2 = 0; ig2 <= ell2; ++ig2) {
			for (ig3 = 0; ig3 <= ell3; ++ig3) {
				pmutm[ig1][ig2][ig3] = 0.0;						/* zero the various allele-frequency counters */
				psel[ig1][ig2][ig3] = 0.0;
				pgtypexp[ig1][ig2][ig3] = 0.0;
				pmutm1[ig1][ig2][ig3] = 0.0;				
				pmutm2[ig1][ig2][ig3] = 0.0;
			}
		}
	}
	
	igen = 0;
	tcount = 0;
	tint = 0.0;
	counter = 0;

	totw = 0.0;
	totwsq = 0.0;

	meanminhet = 0.0; 
	
	totig1 = 0.0;
	totigsq1 = 0.0;
	totwvarig1 = 0.0;

	totig2 = 0.0;
	totigsq2 = 0.0;
	totwvarig2 = 0.0;

	totig3 = 0.0;
	totigsq3 = 0.0;
	totwvarig3 = 0.0;
	
	totw3mom1 = 0.0;
	totw4mom1 = 0.0;
	totw3mom2 = 0.0;
	totw4mom2 = 0.0;
	totw3mom3 = 0.0;
	totw4mom3 = 0.0;
	
	totperform = 0.0;
	totperformsq = 0.0;
	totdevi = 0.0;
	totdevisq = 0.0;
	totabsdevi = 0.0;
	totabsdevisq = 0.0;

	newhigh1 = ell1;
	newlow1 = 0;
	newhigh2 = ell2;
	newlow2 = 0;
	newhigh3 = ell3;
	newlow3 = 0;

	increment = ne / xinc;						/* increment in generations between statistic calculations (set as a fraction of Ne). */



	/* ******************************************************************************************************************************************* */


	/* Iterate the recursion equations to obtain the equilibrium expectations. */

	while (tint < ngens)  										/* iterate until the stopping criterion has been met. */
	{
		igen = igen + 1;



		/* Set the running upper and lower boundaries to the allelic count classes. */

		lowig1 = newlow1 - 1;
		highig1 = newhigh1 + 1;

		if (lowig1 < 0) { lowig1 = 0; }
		if (highig1 > ell1) { highig1 = ell1; }
		if (newlow1 < 0) { newlow1 = 0; }
		if (newhigh1 > ell1) { newhigh1 = ell1; }

		lowig2 = newlow2 - 1;
		highig2 = newhigh2 + 1;

		if (lowig2 < 0) { lowig2 = 0; }
		if (highig2 > ell2) { highig2 = ell2; }
		if (newlow2 < 0) { newlow2 = 0; }
		if (newhigh2 > ell2) { newhigh2 = ell2; }

		lowig3 = newlow3 - 1;
		highig3 = newhigh3 + 1;

		if (lowig3 < 0) { lowig3 = 0; }
		if (highig3 > ell3) { highig3 = ell3; }
		if (newlow3 < 0) { newlow3 = 0; }
		if (newhigh3 > ell3) { newhigh3 = ell3; }


		for (ig1 = lowig1; ig1 <= highig1; ++ig1) {												/* zero the frequencies of the classes */
			for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
				for (ig3 = lowig3; ig3 <= highig3; ++ig3) {
					psel[ig1][ig2][ig3] = 0.0;
				}
			}
		}



		/* Impose selection. */

		meanfit = 0.0;

		for (ig1 = newlow1; ig1 <= newhigh1; ++ig1) {											/* calculate mean fitness */
			for (ig2 = newlow2; ig2 <= newhigh2; ++ig2) {
				for (ig3 = newlow3; ig3 <= newhigh3; ++ig3) {
					meanfit = meanfit + (p0[ig1][ig2][ig3] * relw[ig1][ig2][ig3]); } } }

		for (ig1 = newlow1; ig1 <= newhigh1; ++ig1) {											/* weight the prior genotype frequencies by relative fitness */
			for (ig2 = newlow2; ig2 <= newhigh2; ++ig2) {
				for (ig3 = newlow3; ig3 <= newhigh3; ++ig3) {
					psel[ig1][ig2][ig3] = p0[ig1][ig2][ig3] * relw[ig1][ig2][ig3] / meanfit; } } }




		/* Impose mutation on the genotypic classes. */

		sump = 0.0;

		if (lowig1 == 0) { lowmut1 = 1; }
		else { lowmut1 = lowig1; }
		if (highig1 == ell1) { highmut1 = ell1 - 1; }
		else { highmut1 = highig1; }

		if (lowig2 == 0) { lowmut2 = 1; }
		else { lowmut2 = lowig2; }
		if (highig2 == ell2) { highmut2 = ell2 - 1; }
		else { highmut2 = highig2; }

		if (lowig3 == 0) { lowmut3 = 1; }
		else { lowmut3 = lowig3; }
		if (highig3 == ell3) { highmut3 = ell3 - 1; }
		else { highmut3 = highig3; }


		for (ig1 = lowmut1; ig1 <= highmut1; ++ig1) {											/* first, do the 1-effect classes */
			for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
				for (ig3 = lowig3; ig3 <= highig3; ++ig3) {
					pmutm1[ig1][ig2][ig3] = (mutransup1[ig1] * psel[ig1 - 1][ig2][ig3]) + (mutransdown1[ig1] * psel[ig1 + 1][ig2][ig3]) + (mutrans01[ig1] * psel[ig1][ig2][ig3]); } } }

		if (lowig1 == 0) {
			for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
				for (ig3 = lowig3; ig3 <= highig3; ++ig3) {
					pmutm1[0][ig2][ig3] = (mut01 * psel[0][ig2][ig3]) + (mutdown01 * psel[1][ig2][ig3]); } } }

		if (highig1 == ell1) {
			for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
				for (ig3 = lowig3; ig3 <= highig3; ++ig3) {
					pmutm1[ell1][ig2][ig3] = (mute1 * psel[ell1][ig2][ig3]) + (mutupe1 * psel[ell1 - 1][ig2][ig3]); } } }



		for (ig1 = lowig1; ig1 <= highig1; ++ig1) {												/* next, do the 2-effect classes */
			for (ig2 = lowmut2; ig2 <= highmut2; ++ig2) {
				for (ig3 = lowig3; ig3 <= highig3; ++ig3) {
					pmutm2[ig1][ig2][ig3] = (mutransup2[ig2] * pmutm1[ig1][ig2 - 1][ig3]) + (mutransdown2[ig2] * pmutm1[ig1][ig2 + 1][ig3]) + (mutrans02[ig2] * pmutm1[ig1][ig2][ig3]); } } }

		if (lowig2 == 0) {
			for (ig1 = lowig1; ig1 <= highig1; ++ig1) {
				for (ig3 = lowig3; ig3 <= highig3; ++ig3) {
					pmutm2[ig1][0][ig3] = (mut02 * pmutm1[ig1][0][ig3]) + (mutdown02 * pmutm1[ig1][1][ig3]); } 	} }

		if (highig2 == ell2) {
			for (ig1 = lowig1; ig1 <= highig1; ++ig1) {
				for (ig3 = lowig3; ig3 <= highig3; ++ig3) {
					pmutm2[ig1][ell2][ig3] = (mute2 * pmutm1[ig1][ell2][ig3]) + (mutupe2 * pmutm1[ig1][ell2 - 1][ig3]); } } }



		for (ig1 = lowig1; ig1 <= highig1; ++ig1) {									/* finally, do the 3-effect classes */
			for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
				for (ig3 = lowmut3; ig3 <= highmut3; ++ig3) {
					pmutm[ig1][ig2][ig3] = (mutransup3[ig3] * pmutm2[ig1][ig2][ig3 - 1]) + (mutransdown3[ig3] * pmutm2[ig1][ig2][ig3 + 1]) + (mutrans03[ig3] * pmutm2[ig1][ig2][ig3]);
					sump = sump + pmutm[ig1][ig2][ig3]; 	} } }

		if (lowig3 == 0) {
			for (ig1 = lowig1; ig1 <= highig1; ++ig1) {
				for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
					pmutm[ig1][ig2][0] = (mut03 * pmutm2[ig1][ig2][0]) + (mutdown03 * pmutm2[ig1][ig2][1]);
					sump = sump + pmutm[ig1][ig2][0]; } } }

		if (highig3 == ell3) {
			for (ig1 = lowig1; ig1 <= highig1; ++ig1) {
				for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
					pmutm[ig1][ig2][ell3] = (mute3 * pmutm2[ig1][ig2][ell3]) + (mutupe3 * pmutm2[ig1][ig2][ell3 - 1]);
					sump = sump + pmutm[ig1][ig2][ell3]; } 	} }




		/* Reset the next generation's expected genotype frequencies, and ensure that they sum to 1.0. */

		for (ig1 = lowig1; ig1 <= highig1; ++ig1) {
			for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
				for (ig3 = lowig3; ig3 <= highig3; ++ig3) {
					pgtypexp[ig1][ig2][ig3] = pmutm[ig1][ig2][ig3] / sump;
					p0[ig1][ig2][ig3] = 0.0; } 	} 	}





		/* Sample the population for new genotype frequencies. */

		ntot = ne;
		sump = 0.0;
		
		newlow1 = ell1;
		newhigh1 = 0;
		newlow2 = ell2;
		newhigh2 = 0;
		newlow3 = ell3;
		newhigh3 = 0;
		
		
		for (ig1 = lowig1; ig1 <= highig1; ++ig1) {
			for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
				for (ig3 = lowig3; ig3 <= highig3; ++ig3) {

					if ((pgtypexp[ig1][ig2][ig3] > 0.0) && (ntot > 0))  {
						pp = pgtypexp[ig1][ig2][ig3] / (1.0 - sump);											/* this is the remaining frequency to sample */

						if (pp >= 1.0000000000000) {
							draw = ntot;
							p0[ig1][ig2][ig3] = ((double)draw) / ((double)ne); 	}

						else {
							draw = gsl_ran_binomial_tpe(rand_new, pp, ntot);
							p0[ig1][ig2][ig3] = ((double)draw) / ((double)ne); 	}

						ntot = ntot - draw;
						sump = sump + pgtypexp[ig1][ig2][ig3];

						if (p0[ig1][ig2][ig3] > 0.0) {
							if (ig1 < newlow1) { newlow1 = ig1; }
							if (ig1 > newhigh1) { newhigh1 = ig1; }
							if (ig2 < newlow2) { newlow2 = ig2; }
							if (ig2 > newhigh2) { newhigh2 = ig2; }
							if (ig3 < newlow3) { newlow3 = ig3; }
							if (ig3 > newhigh3) { newhigh3 = ig3; }
							
						}}}}
		}



		

		/* Calculate the summary statistics if the sampling interval is completed. */

		if (igen == increment) {
			igen = 0;
			counter = counter + 1;

			if (counter > burnin) {
				meanw = 0.0;
				meanig1 = 0.0;
				ssqig1 = 0.0;
				meanig2 = 0.0;
				ssqig2 = 0.0;
				meanig3 = 0.0;
				ssqig3 = 0.0;
				w3momig1 = 0.0;
				w3momig2 = 0.0;
				w3momig3 = 0.0;
				w4momig1 = 0.0;
				w4momig2 = 0.0;
				w4momig3 = 0.0;

				for (ig1 = lowig1; ig1 <= highig1; ++ig1) {
					for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
						for (ig3 = lowig3; ig3 <= highig3; ++ig3) {

							sumfreq1[ig1] = sumfreq1[ig1] + p0[ig1][ig2][ig3];							/* get the haplotype frequency distributions */
						sumfreq2[ig2] = sumfreq2[ig2] + p0[ig1][ig2][ig3];
						sumfreq3[ig3] = sumfreq3[ig3] + p0[ig1][ig2][ig3];

						meanw = meanw + (p0[ig1][ig2][ig3] * corwfit[ig1][ig2][ig3]);					/* get the mean fitness */

						meanig1 = meanig1 + (p0[ig1][ig2][ig3] * ((double)ig1));						/* get the first and second moments of the haplotypes */
						ssqig1 = ssqig1 + (p0[ig1][ig2][ig3] * pow(((double)ig1), 2.0));
						
						meanig2 = meanig2 + (p0[ig1][ig2][ig3] * ((double)ig2));
						ssqig2 = ssqig2 + (p0[ig1][ig2][ig3] * pow(((double)ig2), 2.0)); 

						meanig3 = meanig3 + (p0[ig1][ig2][ig3] * ((double)ig3));
						ssqig3 = ssqig3 + (p0[ig1][ig2][ig3] * pow(((double)ig3), 2.0)); 
						}
					}
				}

				for (ig1 = lowig1; ig1 <= highig1; ++ig1) {												/* get the third and fourth moments */
					for (ig2 = lowig2; ig2 <= highig2; ++ig2) {
						for (ig3 = lowig3; ig3 <= highig3; ++ig3) {

							w3momig1 = w3momig1 + (p0[ig1][ig2][ig3] * pow((((double)ig1) - meanig1), 3.0));
							w4momig1 = w4momig1 + (p0[ig1][ig2][ig3] * pow((((double)ig1) - meanig1), 4.0));

							w3momig2 = w3momig2 + (p0[ig1][ig2][ig3] * pow((((double)ig2) - meanig2), 3.0));
							w4momig2 = w4momig2 + (p0[ig1][ig2][ig3] * pow((((double)ig2) - meanig2), 4.0));

							w3momig3 = w3momig3 + (p0[ig1][ig2][ig3] * pow((((double)ig3) - meanig3), 3.0));
							w4momig3 = w4momig3 + (p0[ig1][ig2][ig3] * pow((((double)ig3) - meanig3), 4.0));
						}
					}
				}


				wvarig1 = (ssqig1 - (meanig1 * meanig1));									/* within population variances for numbers of + alleles */
				wvarig2 = (ssqig2 - (meanig2 * meanig2)) ;
				wvarig3 = (ssqig3 - (meanig3 * meanig3));
				
				totw = totw + meanw;
				totwsq = totwsq + pow(meanw, 2.0);
				
				totig1 = totig1 + meanig1;
				totigsq1 = totigsq1 + pow(meanig1, 2.0);
				totwvarig1 = totwvarig1 + wvarig1;
				totw3mom1 = totw3mom1 + w3momig1;
				totw4mom1 = totw4mom1 + w4momig1;

				totig2 = totig2 + meanig2;
				totigsq2 = totigsq2 + pow(meanig2, 2.0);
				totwvarig2 = totwvarig2 + wvarig2;
				totw3mom2 = totw3mom2 + w3momig2;
				totw4mom2 = totw4mom2 + w4momig2;

				totig3 = totig3 + meanig3;
				totigsq3 = totigsq3 + pow(meanig3, 2.0);
				totwvarig3 = totwvarig3 + wvarig3;
				totw3mom3 = totw3mom3 + w3momig3;
				totw4mom3 = totw4mom3 + w4momig3;
				

				tint = tint + 1.0;
				tcount = tcount + 1;

				grandmeanig1 = totig1 / tint;
				varig1 = totwvarig1 / tint;

				grandmeanig2 = totig2 / tint;
				varig2 = totwvarig2 / tint;

				grandmeanig3 = totig3 / tint;
				varig3 = totwvarig3 / tint;

	            perform = ( (meanig1 * sco1 / sco3) + (meanig2 * sco2 / sco3) + meanig3 ) ;
                totperform = totperform + perform;
                totperformsq = totperformsq + pow(perform, 2.0);

				devi = perform - thetaopt;
				totdevi = totdevi + devi;
				totdevisq = totdevisq + pow(devi, 2.0);

				absdevi = abs(perform - thetaopt);
				totabsdevi = totabsdevi + absdevi;
				totabsdevisq = totabsdevisq + pow(absdevi, 2.0);



				totmin0 = 0.0;											/* Calculates the heterozygosity at the class-1 site. */

				for (ig2 = 0; ig2 <= ell2; ++ig2) {						/* Is only relevant if there is one such site, and it is set to be neutral. */
					for (ig3 = 0; ig3 <= ell3; ++ig3) {
						totmin0 = totmin0 + p0[0][ig2][ig3]; }	}

				meanminhet = meanminhet + (2.0 * totmin0 * (1.0 - totmin0));


				
				if (tcount > tintprint) {
					printf("%9ld, %9d, %9.0f, %7.6f, %10.6f, %7.6f, %7.6f, %6d, %4d, %4d, %6d, %4d, %4d, %6d, %4d, %4d, %10.4f, %7.4f, %7.4f, %10.4f\n", (ne*kfac), kfac, tint, (totw / tint), 
					(grandmeanig1 / ((double)ell1)), (grandmeanig2 / ((double)ell2)), (grandmeanig3 / ((double)ell3)),
						newlow1, newhigh1, (newhigh1 - newlow1), newlow2, newhigh2, (newhigh2 - newlow2), newlow3, newhigh3, (newhigh3 - newlow3), pow(varig1, 0.5), pow(varig2, 0.5), pow(varig3, 0.5), perform);
					
					tcount = 0; }

			}
		}								/* ends the summary statistic analysis for this point */

	}									/* ends the loop for generations. */



	/* Calculate the final statistics. */

	for (ig1 = 0; ig1 <= ell1; ++ig1) {							/* these are the steady-state frequencies of the allelic types for the three types of sites */
		meanpfreq1[ig1] = sumfreq1[ig1] / tint; }

	for (ig2 = 0; ig2 <= ell2; ++ig2) {
		meanpfreq2[ig2] = sumfreq2[ig2] / tint; }

	for (ig3 = 0; ig3 <= ell3; ++ig3) {
		meanpfreq3[ig3] = sumfreq3[ig3] / tint; }

	grandmeanw = totw / tint;									/* mean fitness */

	grandmeanig1 = totig1 / tint;								/* mean frequencies for numbers of + alleles at the three site types */
	grandmeanig2 = totig2 / tint;
	grandmeanig3 = totig3 / tint;
	
	varw = (totwsq / tint) - pow(grandmeanw, 2.0);				/* average within-generation variance for fitness, and numbers of + alleles for three site types */
	varig1 = totwvarig1 / tint;
	varig2 = totwvarig2 / tint;
	varig3 = totwvarig3 / tint;

	w3momig1 = totw3mom1 / tint;								/* average within-generation 3rd moments for numbers of + alleles */
	w3momig2 = totw3mom2 / tint;
	w3momig3 = totw3mom3 / tint;

	w4momig1 = totw4mom1 / tint;								/* average within-generation 4th moments for numbers of + alleles */
	w4momig2 = totw4mom2 / tint;
	w4momig3 = totw4mom3 / tint;
	

	meanperform = ((sco1 * grandmeanig1 / sco3) + (sco2 * grandmeanig2 / sco3) + grandmeanig3);		/* mean performance downweighting minor alleles */

	meanminhet = meanminhet / tint;								/* mean heterozygosity observed at the neutral (site type 1) site */

	neuthet = 4.0 * ne * u01 / (1.0 + mutbeta + (ne * u10 * (1.0 + (6.0*mutbeta) + pow(mutbeta, 2.0))));		/* expected neutral heterozygosity if Ne = N */
	
	neestimate = meanminhet * (1 + mutbeta) / u10;				/* Ne estimate derived from neutral heterozygosity */
	
	neestimate = neestimate / (     (4.0*mutbeta) - (meanminhet * (1.0 + (6.0*mutbeta) + pow(mutbeta,2.0)))       );
	
	fracpine = ((double)neestimate) / ((double)ne);				/* estimated Ne based on neutral site vs. ideal expectation */
	
	mean1 = grandmeanig1 / ((double)ell1); 						/* average allele frequencies for three site types */
	mean2 = grandmeanig2 / ((double)ell2);
	mean3 = grandmeanig3 / ((double)ell3);
	
	bne1 = (0.5 / (selco1 / ((double)kfac))) * log(mean1 / ((1.0 - mean1) * mutbeta));		/* Ne expected from Li-Bulmer to account for mean observed allele frequency */
	bne2 = (0.5 / (selco2 / ((double)kfac))) * log (mean2 / ((1.0 - mean2) * mutbeta));
	bne3 = (0.5 / (selco3 / ((double)kfac))) * log(mean3 / ((1.0 - mean3) * mutbeta));
	
	fracbne1 = bne1 / (((double)ne) * kfac);					/* Li-Bulmer Ne relative to ideal value (without interference) */
	fracbne2 = bne2 / (((double) ne) * kfac);
	fracbne3 = bne3 / (((double)ne) * kfac);

	bmean = totig1 / tint;										/* mean, variance, and SD for numbers of + alleles across generations */
    bmnsq = totigsq1 / tint;
    bvar = bmnsq - (bmean * bmean);
    bsd1 = pow(bvar,0.5) / ((double)ell1);
    
    bmean = totig2 / tint;
    bmnsq = totigsq2 / tint;
    bvar = bmnsq - (bmean * bmean);
    bsd2 = pow(bvar,0.5) / ((double)ell2);
   
    bmean = totig3 / tint;
    bmnsq = totigsq3 / tint;
    bvar = bmnsq - (bmean * bmean);
    bsd3 = pow(bvar,0.5) / ((double)ell3);
    
	bmean = meanperform;										/* mean, variance, and SD for for performance, performance - optimum, and absolute deviation */
    bmnsq = totperformsq / tint;
    bvar = bmnsq - (bmean * bmean);
    bsdperform = pow(bvar, 0.5);

	bmean = totdevi / tint;
	meandevi = bmean;
	bmnsq = totdevisq / tint;
	bvar = bmnsq - (bmean * bmean);
	bsddevi = pow(bvar, 0.5);

	bmean = totabsdevi / tint;
	meanabsdevi = bmean;
	bmnsq = totabsdevisq / tint;
	bvar = bmnsq - (bmean * bmean);
	bsdabsdevi = pow(bvar, 0.5);



	fprintf(stream, " %11ld, %11d, %11d, %11d,,"
		" %12.11f, %12.11f, %12.11f, %3.2f,, %12.11f,,"
		" %12.11f, %12.11f ,,  %9d ,,"
		" %17.0f, %11d, %17.0f ,, %12.11f, %12.11f ,,"
		" %12.11f, %12.11f, %12.11f ,, %12.11f, %12.11f, %12.11f ,, %12.11f, %12.11f, %12.11f ,,"
		" %12.11f,, %12.11f, %12.11f ,, %12.0f ,, %10.8f,, %10.8f, %10.8f, %10.8f,,"
		" %10.8f, %10.8f, %10.8f,, %10.8f,, %10.8f, %10.8f,, %10.8f, %10.8f,,"
		" %10.8f, %10.8f, %10.8f,, %10.8f, %10.8f, %10.8f,, %10.8f, %10.8f\n  ",
		(ne*kfac), ell1, ell2, ell3,
		(selco1 / ((double)kfac)), (selco2 / ((double)kfac)), (selco3 / ((double)kfac)), powerexp, thetaopt,
		(u10 / ((double)kfac)), mutbeta, kfac,
		ngens, burnin, (tint*((double)increment)), grandmeanw, varw,
		mean1, mean2, mean3, pseq1, pseq2, pseq3, pow(varig1, 0.5), pow(varig2, 0.5), pow(varig3, 0.5),
		meanperform, meanminhet, neuthet, (((double)neestimate) * ((double)kfac)), fracpine, fracbne1, fracbne2, fracbne3,
		bsd1, bsd2, bsd3, bsdperform, meandevi, bsddevi, meanabsdevi, bsdabsdevi,
		w3momig1, w3momig2, w3momig3, w4momig1, w4momig2, w4momig3, (meandevi/thetaopt), (meanabsdevi/thetaopt)	);

	printf("\n");

	fclose(stream);

}									/* End the set of iterations over all population sizes and mutation rates. */

exit(0);

}





