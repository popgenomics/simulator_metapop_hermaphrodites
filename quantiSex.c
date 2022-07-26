#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_permutation.h>
#define VERSION "03.01.2021"
#define DEPENDENCY "None\n"
#define MAX_NUMBER_OF_INITIAL_NTRL_ALLELES 999	// number of segregating alleles when generating the first parental population
#define RANGE 0.1	// value in [0;1] to modify the current allelic effect between [(1-RANGE) x current_value ; (1+RANGE) * current_value].
#define KRED  "\033[1m\033[31m"
#define KMAG  "\x1B[31m"
#define STOP  "\x1B[0m"

//	gcc quantiSex.c -L/usr/local/lib -lgsl -lgslcblas -lm -Wall -Wextra -Wshadow -Werror -O3 -o quantiSex
//	./quantiSex 100 200 100 10 0.00001 1 0.0001 10 0 0.3 1 0 1 123

typedef struct Deme Deme;
struct Deme{
	int nIndividus;
	long* ntrlLoci;	// contains the neutral alleles for nIndividus (number of individuals within the deme) x 2 (diploids) x nNtlrLoci (number of neutral loci)
	double* quantiLoci; // contains the allelic effects for nIndividus (number of individuals within the deme) x 2 (diploids) x nQuantiLoci (number of quantitative loci)
	double* femaleAllocation; // =sum of the allelic effects over the quantitative loci
	double* maleAllocation; // =(1 - femaleAllocation)
	double fertility; // sum of individual female allocations within the deme
//	int* nOffsprings; // floor(X) + Binom(n=1, p=X-floor(x)) where X = femaleAllocation x fecundity
};

void initializePopulation(gsl_rng *r, Deme* population, const int nDemes, const int maxIndPerDem, const int nNtrlLoci, const int nQuantiLoci, const double z_exp);
void libererMemoirePopulation(Deme* population, const int nDemes);
void setToZero(const int nDemes, int nImmigrants[], int extinctionStatus[], int nProducedSeeds[]);
void initializeNewPopulation(Deme* newPopulation, const int nDemes, const int maxIndPerDem, const int nProducedSeeds[], const int nNtrlLoci, const int nQuantiLoci, const int fecundity, const int extinctionStatus[]);
void weightedSample(gsl_rng* r, const double* liste, const double* weights, double* target, const int sizeOfListe, const int nTrials);
void writeNindividuals(const Deme* population, const int nDemes, const double extinction, const double migrationRate, const int seed);
void checkCommandLine(int argc);
void statisticsPopulations(Deme* population, const int nDemes, const int maxIndPerDem, const int nQuantiLoci, const int fecundity, const double migrationRate, const double extinction, const double dispersal, const int recolonization, const int sexualSystem, const int seed, int time, const double selfingRate, const int colonizationModel, const double global_fst_cm, const double global_fst_coal, const double global_gpst, const double global_D, const double global_Fis, double avg_nc, double avg_Hs, double avg_Htot, double avg_Ho, const double global_Fst_anova, const double global_Fis_anova, const double global_Fit_anova, const double z_exp);
double fstMullon(const int maxIndPerDem, const double extinction, const int recolonization, const double migrationRate);
double fstRousset(const int maxIndPerDem, const double extinction, const int recolonization, const double migrationRate, const int colonizationModel);
void global_stat(Deme* population, const int nDemes, const long nNtrlLoci, double* diff_stats, double* polym_stats);
void QoneQtwoQthree(const Deme* population, const int nDemes, const int nNtrlLoci, const int locus, double* target);
void nc(const Deme* population, const int nDemes, const int nNtrlLoci, const int locus, double* target);
double heteroZ(const double* cont_table_tot);
int valueInArray(const float val, const int sizeArr, const int* array);

void getExtinctionStatus(gsl_rng* r, const int nDemes, const double extinction, int extinctionStatus[]);
void getProducedSeeds(gsl_rng* r, const int nDemes, Deme* population, const int maxIndPerDem,  const int fecundity, int nProducedSeeds[]);
void getImmigrants(gsl_rng* r, const int nProducedSeeds[], const int nDemes, const int extinctionStatus[], const double migrationRate, int nImmigrants[]);


void breeding_step(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const double selfingRate, const double dispersal, const int extinctionStatus[]);
void breeding_step_v2(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const double selfingRate, const double dispersal, const int extinctionStatus[]);

///
void migration_step(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int extinctionStatus[], const int nImmigrants[]); // backward; cloning; doesn't deal with dispersal
void migration_step_v2(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int extinctionStatus[], const int nImmigrants[], const double selfingRate, const double dispersal); // forward; gets mothers and fathers from the metapop; deals with dispersal

void recolonization_step(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int extinctionStatus[], const int recolonization, const int colonizationModel, const double selfingRate, const int maxIndPerDem); // backward; cloning; doesn't deal with dispersal
void recolonization_step_v2(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int extinctionStatus[], const int recolonization, const int colonizationModel, const double dispersal, const double selfingRate, const int maxIndPerDem); // forward; gets mothers and fathers from the metapop; deals with dispersal
///

void replacement(Deme* population, const Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci);
void get_parent_for_dispersal(gsl_rng *r, const int nParents, double deme_target_mother[], double deme_target_father[], int ind_target_mother[], int ind_target_father[], const Deme* population, const int nDemes, const int extinctionStatus[], const double selfingRate, const double dispersal, const int local_deme, const int dispersal_model);
void afficherPopulation(const Deme* population, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int extinctionStatus[], const int generation);
void mutation_step(gsl_rng* r, Deme* population, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const double ntrlMutation, const double quantiMutation);

double z(const int maxIndPerDem, const double selfingRate, const double migrationRate, const double dispersal, const double extinction, const int recolonization);

int main(int argc, char *argv[]){

	checkCommandLine(argc); // stop the code if the number of arguments doesn't fit with the expected one

	int i = 0;
	int j = 0;

	// Get Parameters from comamnd line
	const int nDemes = atoi(argv[1]); // number of demes
	const int maxIndPerDem = atoi(argv[2]);	// carrying capacity per deme
	const int nGeneration = atoi(argv[3]);	// number of generations to simulate
//	const int nGenerationUnisex = atoi(argv[4]);	// generations at which unisexuals appear in the metapopulation
	const int nNtrlLoci = atoi(argv[4]);	// number of neutral loci
	const double ntrlMutation = atof(argv[5]);	// mutation rate of the ntrl loci
	const int nQuantiLoci = atoi(argv[6]);	// number of quantitative loci
	const double quantiMutation = atof(argv[7]);	// mutation rate of the quantative loci
	const int fecundity = atoi(argv[8]);	// max number of offspring when femAlloc=100%
	const double migrationRate = atof(argv[9]);	// immigration rate
	const double dispersal = atof(argv[10]);	// pollen dispersal rate 
	const double extinction = atof(argv[11]);	// extinction rate
	const int recolonization = atoi(argv[12]);	// number of recolonizing individuals
	const int colonizationModel = atoi(argv[13]);	// 0 = migrant pool model 1 = propagule pool model
	const int sexualSystem = atoi(argv[14]);    // 0 = only hermaphrodites; 1 = XY system; 2 = ZW system
//	const double sexAvantage = atof(argv[16]); // avantage confered by the Y or Z chromosome over hermaphrodites
	const double selfingRate = atof(argv[15]); // probability to have an ovule being fertilized by sperm from the same individual
	const int lapse_stats = atoi(argv[16]); // print statistics every "lapse_stats" generations
	const int seed = atoi(argv[17]);

	// expected z
	double z_exp = z(maxIndPerDem, selfingRate,  migrationRate, dispersal, extinction, recolonization);
	
	// Random generator
        const gsl_rng_type *T;
        gsl_rng *r;
        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, seed);

	// array of statistics summaryzing the genetic differentiation
	double* diff_stats = NULL;
	double* polym_stats = NULL;
	diff_stats = malloc(8 * sizeof(double));
	polym_stats = malloc(4 * sizeof(double));
	if( diff_stats == NULL || polym_stats == NULL){
		exit(0);
	}	

	// Initializing the metapopulation
	Deme* population = NULL;
	population = malloc(nDemes * sizeof(Deme));
	if(population == NULL){
		exit(0);
	}

	initializePopulation(r, population, nDemes, maxIndPerDem, nNtrlLoci, nQuantiLoci, z_exp);

	// fst, gst, gpst, Jost D
	double global_fst_cm = 0.0;
	double global_fst_coal = 0.0;
	double global_gpst = 0.0;
	double global_D = 0.0;
	double global_Fis = 0.0;

	double global_Fst_anova = 0.0;
	double global_Fis_anova = 0.0;
	double global_Fit_anova = 0.0;
	
	global_Fst_anova+=0;
	global_Fis_anova+=0;
	global_Fit_anova+=0;

	
	// polymorphisms
	double avg_nc = 0.0;
	double avg_Hs = 0.0;
	double avg_Htot = 0.0;
	double avg_Ho = 0.0;

	// Evolution of the metapopulation
	for(i=0; i<=nGeneration; i++){	// start of the loop 'i' over the generations
		int* nImmigrants = NULL; // stock the number of received immigrants per deme
		int* extinctionStatus = NULL;	// per deme: 0 = non-extincted; 1 = extincted
		int* nProducedSeeds = NULL; // stock the numer of produced seeds per deme
		nImmigrants = malloc(nDemes * sizeof(int));
		extinctionStatus = malloc(nDemes * sizeof(int));
		nProducedSeeds = malloc(nDemes * sizeof(int));
		
		//afficherPopulation(population, nDemes, nNtrlLoci, nQuantiLoci, extinctionStatus, i);
		
		if(nImmigrants == NULL || extinctionStatus == NULL || nProducedSeeds == NULL){
			exit(0);
		}
		setToZero(nDemes, nImmigrants, extinctionStatus, nProducedSeeds);	// initialize vectors used to configure the new-population

		// get the extinction status
		getExtinctionStatus(r, nDemes, extinction, extinctionStatus);
		
		// get the number of produced seeds
		getProducedSeeds(r, nDemes, population, maxIndPerDem, fecundity, nProducedSeeds);
		
		// get the number of immigrants
		getImmigrants(r, nProducedSeeds, nDemes, extinctionStatus, migrationRate, nImmigrants);
		
		// new population = seeds (not germinated yet)
		Deme* newPopulation = NULL;
		newPopulation = malloc(nDemes * sizeof(Deme));
		if(newPopulation == NULL){
			exit(0);
		}
		initializeNewPopulation(newPopulation, nDemes, maxIndPerDem, nProducedSeeds, nNtrlLoci, nQuantiLoci, fecundity, extinctionStatus);

//		afficherPopulation(population, nDemes, nNtrlLoci, nQuantiLoci, extinctionStatus, i);
		// reproduction
//		breeding_step(r, population, newPopulation, nDemes, nNtrlLoci, nQuantiLoci, selfingRate, dispersal, extinctionStatus); // fast: conditionnal pollen dispersal
		breeding_step_v2(r, population, newPopulation, nDemes, nNtrlLoci, nQuantiLoci, selfingRate, dispersal, extinctionStatus); // slow: large fitness vector
		
		// migration
		migration_step(r, population, newPopulation, nDemes, nNtrlLoci, nQuantiLoci, extinctionStatus, nImmigrants); // backward 
//		migration_step_v2(r, population, newPopulation, nDemes, nNtrlLoci, nQuantiLoci, extinctionStatus, nImmigrants, selfingRate, dispersal); // forward 
		
		// recolonization
		
		recolonization_step(r, population, newPopulation, nDemes, nNtrlLoci, nQuantiLoci, extinctionStatus, recolonization, colonizationModel, selfingRate, maxIndPerDem); // cloning
//		recolonization_step_v2(r, population, newPopulation, nDemes, nNtrlLoci, nQuantiLoci, extinctionStatus, recolonization, colonizationModel, dispersal, selfingRate, maxIndPerDem); // not cloning

		// mutation
		mutation_step(r, newPopulation, nDemes, nNtrlLoci, nQuantiLoci, ntrlMutation, quantiMutation);
		
		// compute some stats
		if( i%lapse_stats==0 && i>0 ){
			for(j=0; j<8; ++j){
				diff_stats[j] = 0.0;
			}
			
			for(j=0; j<4; ++j){
				polym_stats[j] = 0.0;
			}
			
			global_stat(newPopulation, nDemes, nNtrlLoci, diff_stats, polym_stats);
			global_fst_cm = diff_stats[0];
			global_fst_coal = diff_stats[1];
			global_gpst = diff_stats[2];
			global_D = diff_stats[3];
			global_Fis = diff_stats[4];

			global_Fst_anova = diff_stats[5];
			global_Fis_anova = diff_stats[6];
			global_Fit_anova = diff_stats[7];

			avg_nc = polym_stats[0];
			avg_Hs = polym_stats[1];
			avg_Htot = polym_stats[2];
			avg_Ho = polym_stats[3];
			
			statisticsPopulations(newPopulation, nDemes, maxIndPerDem, nQuantiLoci, fecundity, migrationRate, extinction, dispersal, recolonization, sexualSystem, seed, i, selfingRate, colonizationModel, global_fst_cm, global_fst_coal, global_gpst, global_D, global_Fis, avg_nc, avg_Hs, avg_Htot, avg_Ho, global_Fst_anova, global_Fis_anova, global_Fit_anova, z_exp);
		}

		//statisticsPopulations(newPopulation, nDemes, maxIndPerDem, nQuantiLoci, fecundity, migration, extinction, recolonization, sexualSystem, sexAvantage, seed, i, selfingRate, colonizationModel, fst_mean);
		
		//writeNindividuals(newPopulation, nDemes, extinction, migration, seed); // to un-comment only if we want #individuals and femaleAllocation per deme per generation

		// remplacer population par newPopulation
		// Nettoyer memoire
		libererMemoirePopulation(population, nDemes);
		free(population);
		population = malloc(nDemes * sizeof(Deme));

		if(population == NULL){
			exit(0);
		}

		replacement(population, newPopulation, nDemes, nNtrlLoci, nQuantiLoci); // replace the parents (poplation) by the offspring (newPopulation)
		
		free(nImmigrants);
		free(extinctionStatus);
		free(nProducedSeeds);
		libererMemoirePopulation(newPopulation, nDemes);
		free(newPopulation);
		
	}	// end of the loop 'i' over the generations
	
	libererMemoirePopulation(population, nDemes);
	
	free(population);
	
	free(polym_stats);	
	free(diff_stats);
	
	return(0);
}


void initializePopulation(gsl_rng* r, Deme* population, const int nDemes, const int maxIndPerDem, const int nNtrlLoci, const int nQuantiLoci, const double z_exp){
	int valuesNtrlAlleles = MAX_NUMBER_OF_INITIAL_NTRL_ALLELES;
	double minQuanti = 0.0;
	double maxQuanti = 0.0;
	maxQuanti = 1.0/2/nQuantiLoci;
//	maxQuanti = 0.5/2/nQuantiLoci;	// uncomment to fix sex allocation to 0.5
	int i = 0;
	for(i=0; i<nDemes; i++){	// loop along the demes
		population[i].nIndividus = maxIndPerDem;
		population[i].ntrlLoci = malloc(2 * maxIndPerDem * nNtrlLoci * sizeof(long));
		population[i].quantiLoci = malloc(2 * maxIndPerDem * nQuantiLoci * sizeof(double));
		population[i].femaleAllocation = malloc(maxIndPerDem * sizeof(double));
		population[i].maleAllocation = malloc(maxIndPerDem * sizeof(double));
		population[i].fertility = 0.0;

		if(population[i].ntrlLoci == NULL || population[i].quantiLoci == NULL || population[i].femaleAllocation == NULL || population[i].maleAllocation == NULL){
			exit(0);
		}

		int cnt = -1;
		int j = 0;
		for(j=0; j<maxIndPerDem; j++){	// loop along the individuals
			cnt += 1;
			population[i].femaleAllocation[j] = 0.0;
			int k = 0;
			for(k=0; k<(2*nNtrlLoci); k++){	// loop along the {2: diploid} x {nNtrlLoci: number of neutral loci} positions
				population[i].ntrlLoci[j*2*nNtrlLoci+k] = gsl_rng_uniform_int(r, valuesNtrlAlleles) + 1;
			}
			for(k=0; k<(2*nQuantiLoci); k++){	// loop along the {2: diploid} x {nQuantiLoci: number of quantitative loci} positions
				population[i].quantiLoci[j*2*nQuantiLoci+k] = gsl_ran_flat(r, minQuanti, maxQuanti);
				population[i].quantiLoci[j*2*nQuantiLoci+k] = z_exp/2/nQuantiLoci; // start from equilibrium
//				population[i].quantiLoci[j*2*nQuantiLoci+k] = maxQuanti; // uncomment to fix sex allocation to 0.5
				population[i].femaleAllocation[j] += population[i].quantiLoci[j*2*nQuantiLoci+k];
//				population[i].maleAllocation[j] -= population[i].quantiLoci[j*2*nQuantiLoci+k];
			}
			
			population[i].maleAllocation[j] = 1-population[i].femaleAllocation[j];
			population[i].fertility += population[i].femaleAllocation[j];

//            		population[i].nOffsprings[j] = floor(fecundity * population[i].femaleAllocation[j]) + gsl_ran_binomial(r, (fecundity * population[i].femaleAllocation[j]) - floor(fecundity * population[i].femaleAllocation[j]), 1);	// nOffs = floor(fecundity x femaleAllocation) + 1 according to a random Binomial integer

        }	// end of loop along the individuals
    }	// end of loop along the demes
}

void setToZero(const int nDemes, int nImmigrants[], int extinctionStatus[], int nProducedSeeds[]){
	// function to set all to demes to zero for the #of immigrants, the extinction status and the #of produced seeds.
	// before configMetapop()
	int i = 0;
	for(i=0; i<nDemes; i++){
		nImmigrants[i] = 0;
		extinctionStatus[i] = 0;
		nProducedSeeds[i] = 0;
	}
}

void getExtinctionStatus(gsl_rng* r, const int nDemes, const double extinction, int extinctionStatus[]){
	int i = 0;
	const unsigned int binomTrials = 1;
	for(i=0; i<nDemes; i++){
		extinctionStatus[i] = gsl_ran_binomial(r, extinction, binomTrials);
	}
}

void getProducedSeeds(gsl_rng* r, const int nDemes, Deme* population, const int maxIndPerDem,  const int fecundity, int nProducedSeeds[]){
	int i = 0;
	int j = 0;
	for(i=0; i<nDemes; i++){
		nProducedSeeds[i] = 0;
		for(j=0; j<population[i].nIndividus; j++){
			nProducedSeeds[i] += floor(fecundity * population[i].femaleAllocation[j]) + gsl_ran_binomial(r, (fecundity * population[i].femaleAllocation[j]) - floor(fecundity * population[i].femaleAllocation[j]), 1); // nOffs = floor(fecundity x femaleAllocation) + 1 according to a random Binomial integer
		}
	
		if(nProducedSeeds[i] > maxIndPerDem){
			nProducedSeeds[i] = maxIndPerDem;
		}
	}
}


void getImmigrants(gsl_rng* r, const int nProducedSeeds[], const int nDemes, const int extinctionStatus[], const double migrationRate, int nImmigrants[]){
	int deme_tmp = 0;
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==1){ // if extinction
			nImmigrants[deme_tmp] = 0;
		}else{
			nImmigrants[deme_tmp] = gsl_ran_poisson(r, migrationRate);
			if(nImmigrants[deme_tmp] > nProducedSeeds[deme_tmp]){
				nImmigrants[deme_tmp] = nProducedSeeds[deme_tmp];
			}
		}
	}
}


void initializeNewPopulation(Deme* newPopulation, const int nDemes, const int maxIndPerDem, const int nProducedSeeds[], const int nNtrlLoci, const int nQuantiLoci, const int fecundity, const int extinctionStatus[]){
	// function to initialize the new population: allocate memory and set values to 0
	int i = 0;
	for(i=0; i<nDemes; i++){ // start of the loop along demes
		int taille = 0;
		if(extinctionStatus[i]==1){
//			taille = recolonization;
			taille = maxIndPerDem;
		}else{
			taille = nProducedSeeds[i];	// size of the deme = nProducedSeeds within the parental population
			if(taille < fecundity){	// if not enough seeds are produced ==> deme is considered as extincted
				taille = fecundity;
			}
			if(taille > maxIndPerDem){	// if too many individuals have to be present in the deme ==> cutoff to the maxIndPerDem (=carrying capacity)
				taille = maxIndPerDem;
			}
		}

		newPopulation[i].nIndividus = taille;
		newPopulation[i].fertility = 0.0;
		
		newPopulation[i].ntrlLoci = NULL;
		newPopulation[i].quantiLoci = NULL;
		newPopulation[i].femaleAllocation = NULL;
		newPopulation[i].maleAllocation = NULL;

		newPopulation[i].ntrlLoci = malloc(2 * taille * nNtrlLoci * sizeof(long));
		newPopulation[i].quantiLoci = malloc(2 * taille * nQuantiLoci * sizeof(long));
		newPopulation[i].femaleAllocation = malloc(taille * sizeof(double));
		newPopulation[i].maleAllocation = malloc(taille * sizeof(double));

		if(newPopulation[i].ntrlLoci == NULL || newPopulation[i].quantiLoci == NULL || newPopulation[i].femaleAllocation == NULL || newPopulation[i].maleAllocation == NULL){
			exit(0);
		}

		int j = 0;
		for(j=0; j<taille; j++){ // start the loop along individuals
			newPopulation[i].femaleAllocation[j] = 0.0;
			newPopulation[i].maleAllocation[j] = 0.0;

			int k = 0;
			for(k=0; k<(2*nNtrlLoci); k++){ // loop along the {2: diploid} x {nNtrlLoci: number of neutral loci} positions
				newPopulation[i].ntrlLoci[j*2*nNtrlLoci + k] = 0;
			}
			for(k=0; k<(2*nQuantiLoci); k++){       // loop along the {2: diploid} x {nQuantiLoci: number of quantitative loci} positions
				newPopulation[i].quantiLoci[j*2*nQuantiLoci + k] = 0;
			}
		}	// end of the loop along individuals
	} // end of the loop along demes
}

void breeding_step(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const double selfingRate, const double dispersal, const int extinctionStatus[]){
	// get fecundity about all demes
	int deme_tmp=0;
	int nNonExtinctedDemes=0;
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==0){
			nNonExtinctedDemes+=1;
		}
	}
	
	double* list_nonExtincted_demes_IDs = NULL;
	list_nonExtincted_demes_IDs = malloc(nNonExtinctedDemes*sizeof(double));
	double* list_nonExtincted_demes_fertility = NULL;
	list_nonExtincted_demes_fertility = malloc(nNonExtinctedDemes*sizeof(double));
	
	if(list_nonExtincted_demes_IDs==NULL || list_nonExtincted_demes_fertility==NULL){
		exit(0);
	}
	
	int cnt=0;
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==0){
			list_nonExtincted_demes_IDs[cnt] = deme_tmp;
			list_nonExtincted_demes_fertility[cnt] = 1-population[deme_tmp].fertility;
			cnt+=1;
		}
	}
	
	// loop over demes
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==0){
			// adjust the male probabilities of being father 
			// and the female probabilities of being a mother 
			double* maleFitness = NULL;
			maleFitness = malloc(population[deme_tmp].nIndividus * sizeof(double));
			
			if(maleFitness == NULL){
				exit(0);
			}
			
			double* femaleFitness = NULL;
			femaleFitness = malloc(population[deme_tmp].nIndividus * sizeof(double));
			if(femaleFitness == NULL){
				exit(0);
			}
			
			double* indIDs = NULL;
			indIDs = malloc(population[deme_tmp].nIndividus * sizeof(double));
			if(indIDs == NULL){
				exit(0);
			}
			
			// get the male and female fitness, and IDs
			int ind_tmp = 0;
			for(ind_tmp=0; ind_tmp<population[deme_tmp].nIndividus; ind_tmp++){
				maleFitness[ind_tmp] = newPopulation[deme_tmp].maleAllocation[ind_tmp];
				femaleFitness[ind_tmp] = newPopulation[deme_tmp].femaleAllocation[ind_tmp];
				indIDs[ind_tmp] = ind_tmp;
			}
			
			
			// get the parents
			double* mothers_ind = NULL;	// array containing the mothers ID. Ex: {2, 19, 3, 3, 1} means that the first baby has individual#2 as mother, second baby has individual#19 as mother, babies 3 and 4 have individual#3
			double* fathers_ind = NULL;	// array containing the fathers ID.
			
			double* mothers_deme = NULL;	// array containing the mothers ID. Ex: {2, 19, 3, 3, 1} means that the first baby has individual#2 as mother, second baby has individual#19 as mother, babies 3 and 4 have individual#3
			double* fathers_deme = NULL;	// array containing the fathers ID.
			
			
			mothers_ind = malloc(newPopulation[deme_tmp].nIndividus * sizeof(double));
			fathers_ind = malloc(newPopulation[deme_tmp].nIndividus * sizeof(double));
			mothers_deme = malloc(newPopulation[deme_tmp].nIndividus * sizeof(double));
			fathers_deme = malloc(newPopulation[deme_tmp].nIndividus * sizeof(double));
			
			if(mothers_ind == NULL || fathers_ind == NULL || mothers_deme == NULL || fathers_deme == NULL ){
				exit(0);
			}
			
			// 1) get the mothers
			weightedSample(r, indIDs, femaleFitness, mothers_ind, population[deme_tmp].nIndividus, newPopulation[deme_tmp].nIndividus);
			
			// 2) get the fathers
			weightedSample(r, indIDs,  maleFitness,  fathers_ind, population[deme_tmp].nIndividus, newPopulation[deme_tmp].nIndividus);
			
			free(femaleFitness);
			free(maleFitness);
			free(indIDs);
			
			// 3) get the demes
			for(ind_tmp=0; ind_tmp<newPopulation[deme_tmp].nIndividus; ind_tmp++){
				mothers_deme[ind_tmp] = deme_tmp;
				fathers_deme[ind_tmp] = deme_tmp;
			}
			
			// 4) deals with pollen dispersal
			for(ind_tmp=0; ind_tmp<newPopulation[deme_tmp].nIndividus; ind_tmp++){
				int test_dispersal=0;
				test_dispersal = gsl_ran_binomial(r, dispersal, 1);
				if(test_dispersal==1){
					// get the deme sending the pollen
					double* dest_deme = NULL;
					dest_deme = malloc(1*sizeof(double));
					if(dest_deme==NULL){
						exit(0);
					}
					weightedSample(r, list_nonExtincted_demes_IDs, list_nonExtincted_demes_fertility, dest_deme, nNonExtinctedDemes, 1);
					
					fathers_deme[ind_tmp]=dest_deme[0];
					free(dest_deme);
					
					// get the male
					double* out_pollen_ID = NULL;
					out_pollen_ID = malloc(population[ (int)fathers_deme[ind_tmp] ].nIndividus * sizeof(double));
					if(out_pollen_ID == NULL){
						exit(0);
					}
					
					double* out_pollen_fitness = NULL;
					out_pollen_fitness = malloc(population[ (int)fathers_deme[ind_tmp] ].nIndividus * sizeof(double));
					if(out_pollen_fitness == NULL){
						exit(0);
					}
					
					int i=0;
					for(i=0; i<population[(int)fathers_deme[ind_tmp] ].nIndividus; i++){
						out_pollen_ID[i] = i;
						out_pollen_fitness[i] = population[(int)fathers_deme[ind_tmp] ].maleAllocation[i];
					}
					
					double* dest_ind = NULL;
					dest_ind = malloc(1*sizeof(double));
					if(dest_ind == NULL){
						exit(0);
					}
					weightedSample(r, out_pollen_ID, out_pollen_fitness, dest_ind, population[(int)fathers_deme[ind_tmp] ].nIndividus, 1);
					
					fathers_ind[ind_tmp]=dest_ind[0];
					free(out_pollen_ID);
					free(out_pollen_fitness);
					free(dest_ind);
				}
			}
				
			// 5) adjust for selfing
			for(ind_tmp=0; ind_tmp<newPopulation[deme_tmp].nIndividus; ind_tmp++){
				int autofec = 0;
				autofec = gsl_ran_binomial(r, selfingRate, 1);
				if(autofec==1){
					fathers_deme[ind_tmp] = deme_tmp;
					fathers_ind[ind_tmp] = mothers_ind[ind_tmp];
				}
			}
			
			// 4) produce the seeds
			int seed_tmp = 0;
			for(seed_tmp=0; seed_tmp<newPopulation[deme_tmp].nIndividus; seed_tmp++){
				// neutral loci
				int ntrlLoci_tmp = 0;
				for(ntrlLoci_tmp=0; ntrlLoci_tmp<nNtrlLoci; ntrlLoci_tmp++){ // loop along the neutral loci
					// mother
					int pos_seed = 0; // position within the receiving deme
					int pos_parent = 0; // position within the donor deme
					pos_seed = 2*nNtrlLoci*seed_tmp + 2*ntrlLoci_tmp + 0;
					pos_parent = 2*nNtrlLoci*(int)mothers_ind[seed_tmp] + 2*ntrlLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].ntrlLoci[pos_seed] = population[(int)mothers_deme[seed_tmp]].ntrlLoci[pos_parent];
					
					// father
					pos_seed += 1;
					pos_parent = 2*nNtrlLoci*(int)fathers_ind[seed_tmp] + 2*ntrlLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].ntrlLoci[pos_seed] = population[(int)fathers_deme[seed_tmp]].ntrlLoci[pos_parent];
					
				} // end of loop for neutral loci
				
				// quantitative loci
				int quantiLoci_tmp = 0;
				for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){ // loop along the quantitative loci
					// mother
					int pos_seed = 0; // position within the receiving deme
					int pos_parent = 0; // position within the donor deme
					pos_seed = 2*nQuantiLoci*seed_tmp + 2*quantiLoci_tmp + 0;
					pos_parent = 2*nQuantiLoci*(int)mothers_ind[seed_tmp] + 2*quantiLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].quantiLoci[pos_seed] = population[deme_tmp].quantiLoci[pos_parent];
					
					// father
					pos_seed += 1;
					pos_parent = 2*nQuantiLoci*(int)fathers_ind[seed_tmp] + 2*quantiLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].quantiLoci[pos_seed] = population[(int)fathers_deme[seed_tmp]].quantiLoci[pos_parent];
				} // end of loop for quantitative loci
			} // end of loop to produce seeds
			
			
			// sex allocation
			newPopulation[deme_tmp].fertility = 0.0;
			for(seed_tmp=0; seed_tmp<newPopulation[deme_tmp].nIndividus; seed_tmp++){
				newPopulation[deme_tmp].femaleAllocation[seed_tmp] = 0.0;
				
				int quantiLoci_tmp = 0;
				for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){ // loop along the quantitative loci
					newPopulation[deme_tmp].femaleAllocation[seed_tmp] += newPopulation[deme_tmp].quantiLoci[2*nQuantiLoci*seed_tmp + 2*quantiLoci_tmp + 0];
					newPopulation[deme_tmp].femaleAllocation[seed_tmp] += newPopulation[deme_tmp].quantiLoci[2*nQuantiLoci*seed_tmp + 2*quantiLoci_tmp + 1];
				}
				newPopulation[deme_tmp].maleAllocation[seed_tmp] = 1 - newPopulation[deme_tmp].femaleAllocation[seed_tmp];
				
				newPopulation[deme_tmp].fertility += newPopulation[deme_tmp].femaleAllocation[seed_tmp];
			}
			
			free(mothers_ind);
			free(fathers_ind);
			free(mothers_deme);
			free(fathers_deme);
		} // end of condition based on the extinction status
	} // end of loop over demes
	
	// end of function
	free(list_nonExtincted_demes_IDs);
	free(list_nonExtincted_demes_fertility);
}



void breeding_step_v2(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const double selfingRate, const double dispersal, const int extinctionStatus[]){
	// get the total number of possible fathers
	int deme_tmp = 0;
	int nIndividusTotal = 0;
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==0){ // if the deme is not extincted
			nIndividusTotal += population[deme_tmp].nIndividus;
		}
	}
	
	// get informations about all possible fathers
	int* deme_ind = NULL;
	int* IDind_demes = NULL;
	double* IDind_metapop = NULL;
	double* maleAlloc_metapop = NULL;
	double* femaleAlloc_metapop = NULL;
	
	deme_ind = malloc(nIndividusTotal * sizeof(int)); // [deme 0; deme 0; deme 0; deme 1; deme 1; deme 1]
	IDind_demes = malloc(nIndividusTotal * sizeof(int)); // [0; 1; 2; 0; 1; 2]
	IDind_metapop = malloc(nIndividusTotal * sizeof(double)); // [0; 1; 2; 3; 4; 5]
	maleAlloc_metapop = malloc(nIndividusTotal * sizeof(double)); // [0.6; 0.7; 0.5; 0.5; 0.6; 0.9]
	femaleAlloc_metapop = malloc(nIndividusTotal * sizeof(double)); // [0.4; 0.3; 0.5; 0.5; 0.4; 0.1]

	if(deme_ind == NULL || IDind_demes == NULL || IDind_metapop == NULL || maleAlloc_metapop == NULL || femaleAlloc_metapop == NULL){
		exit(0);
	}
	
	int cnt = 0;
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==0){
			int ind_tmp = 0;
			for(ind_tmp=0; ind_tmp<population[deme_tmp].nIndividus; ind_tmp++){
				deme_ind[cnt] = deme_tmp; // deme's id for each male
				IDind_demes[cnt] = ind_tmp; // male's id within its deme
				IDind_metapop[cnt] = (double)cnt; // male's id within the metapopulation
				maleAlloc_metapop[cnt] = population[deme_tmp].maleAllocation[ind_tmp]; // male's allocation
				femaleAlloc_metapop[cnt] = population[deme_tmp].femaleAllocation[ind_tmp]; // female's allocation
				cnt += 1;
			}
		}
	}
	
	// loop over demes
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==0){
			// adjust the male probabilities of being father 
			// and the female probabilities of being a mother 
			double* maleFitness = NULL;
			maleFitness = malloc(nIndividusTotal * sizeof(double));
			
			if(maleFitness == NULL){
				exit(0);
			}
			
			// get the fitness for males and females
			int ind_tmp = 0;
			for(ind_tmp=0; ind_tmp<nIndividusTotal; ind_tmp++){
				if( extinctionStatus[deme_ind[ind_tmp]]==0 ){ // don't deal with extincted demes
					if(deme_ind[ind_tmp]==deme_tmp){ // if the male is from the current deme (deme_tmp)
						maleFitness[ind_tmp] = maleAlloc_metapop[ind_tmp]*(1-dispersal);
					}else{ // if the male is from a deme different of deme_tmp
						maleFitness[ind_tmp] = maleAlloc_metapop[ind_tmp]*dispersal / (nDemes-1);
					}
				}
			}
			
			double* femaleFitness = NULL;
			double* femaleIDs = NULL;
			femaleFitness = malloc(population[deme_tmp].nIndividus * sizeof(double));
			femaleIDs = malloc(population[deme_tmp].nIndividus * sizeof(double));
			if(maleFitness == NULL){
				exit(0);
			}
			
			for(ind_tmp=0; ind_tmp<population[deme_tmp].nIndividus; ind_tmp++){
				femaleFitness[ind_tmp] = population[deme_tmp].femaleAllocation[ind_tmp];
				femaleIDs[ind_tmp] = ind_tmp;
			}
			
			// get the parents
			double* mothers = NULL;	// array containing the mothers ID. Ex: {2, 19, 3, 3, 1} means that the first baby has individual#2 as mother, second baby has individual#19 as mother, babies 3 and 4 have individual#3
			double* fathers = NULL;	// array containing the fathers ID.
			
			mothers = malloc(newPopulation[deme_tmp].nIndividus * sizeof(double));
			fathers = malloc(newPopulation[deme_tmp].nIndividus * sizeof(double));
			
			if(mothers == NULL || fathers == NULL){
				exit(0);
			}
			
			// 1) get the mothers
			weightedSample(r, femaleIDs, femaleFitness, mothers, population[deme_tmp].nIndividus, newPopulation[deme_tmp].nIndividus);
			
			// 2) get the fathers
//			printf("%d\t%d\n", nIndividusTotal, newPopulation[deme_tmp].nIndividus);
			weightedSample(r, IDind_metapop, maleFitness,   fathers, nIndividusTotal, newPopulation[deme_tmp].nIndividus);
			
			// 3) adjust for selfing
			double* fathers_demes = NULL;
			fathers_demes = malloc(newPopulation[deme_tmp].nIndividus*sizeof(double));
			if(fathers_demes==NULL){
				exit(0);
			}
			
			for(ind_tmp=0; ind_tmp<newPopulation[deme_tmp].nIndividus; ind_tmp++){
				int autofec = 0;
				autofec = gsl_ran_binomial(r, selfingRate, 1);
				if(autofec==1){
					fathers_demes[ind_tmp] = deme_tmp;
					fathers[ind_tmp] = mothers[ind_tmp];
				}else{
					fathers_demes[ind_tmp] = deme_ind[(int)fathers[ind_tmp]];
				}
			}
			
			// 4) produce the seeds
			int seed_tmp = 0;
			for(seed_tmp=0; seed_tmp<newPopulation[deme_tmp].nIndividus; seed_tmp++){
				// neutral loci
				int ntrlLoci_tmp = 0;
				for(ntrlLoci_tmp=0; ntrlLoci_tmp<nNtrlLoci; ntrlLoci_tmp++){ // loop along the neutral loci
					// mother
					int pos_seed = 0; // position within the receiving deme
					int pos_parent = 0; // position within the donor deme
					pos_seed = 2*nNtrlLoci*seed_tmp + 2*ntrlLoci_tmp + 0;
					pos_parent = 2*nNtrlLoci*mothers[seed_tmp] + 2*ntrlLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].ntrlLoci[pos_seed] = population[deme_tmp].ntrlLoci[pos_parent];
					
					// father
					pos_seed += 1;
					pos_parent = 2*nNtrlLoci*IDind_demes[(int)fathers[seed_tmp]] + 2*ntrlLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].ntrlLoci[pos_seed] = population[(int)fathers_demes[seed_tmp]].ntrlLoci[pos_parent];
				} // end of loop for neutral loci
				
				// quantitative loci
				int quantiLoci_tmp = 0;
				for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){ // loop along the quantitative loci
					// mother
					int pos_seed = 0; // position within the receiving deme
					int pos_parent = 0; // position within the donor deme
					pos_seed = 2*nQuantiLoci*seed_tmp + 2*quantiLoci_tmp + 0;
					pos_parent = 2*nQuantiLoci*mothers[seed_tmp] + 2*quantiLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].quantiLoci[pos_seed] = population[deme_tmp].quantiLoci[pos_parent];
					
					// father
					pos_seed += 1;
					pos_parent = 2*nQuantiLoci*IDind_demes[(int)fathers[seed_tmp]] + 2*quantiLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].quantiLoci[pos_seed] = population[(int)fathers_demes[seed_tmp]].quantiLoci[pos_parent];
				} // end of loop for quantitative loci
			} // end of loop to produce seeds
			
			
			// sex allocation
			for(seed_tmp=0; seed_tmp<newPopulation[deme_tmp].nIndividus; seed_tmp++){
				int quantiLoci_tmp = 0;
				newPopulation[deme_tmp].femaleAllocation[seed_tmp] = 0.0;
				for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){ // loop along the quantitative loci
					newPopulation[deme_tmp].femaleAllocation[seed_tmp] += newPopulation[deme_tmp].quantiLoci[2*nQuantiLoci*seed_tmp + 2*quantiLoci_tmp + 0];
					newPopulation[deme_tmp].femaleAllocation[seed_tmp] += newPopulation[deme_tmp].quantiLoci[2*nQuantiLoci*seed_tmp + 2*quantiLoci_tmp + 1];
				}
				newPopulation[deme_tmp].maleAllocation[seed_tmp] = 1 - newPopulation[deme_tmp].femaleAllocation[seed_tmp];
				
				newPopulation[deme_tmp].fertility += newPopulation[deme_tmp].femaleAllocation[seed_tmp];
			}
			
			free(mothers);
			free(fathers);
			free(maleFitness);
			free(femaleFitness);
			free(femaleIDs);
			free(fathers_demes);
		} // end of condition based on the extinction status
	} // end of loop over demes
	
	// end of function
	free(deme_ind); // [deme 0; deme 0; deme 0; deme 1; deme 1; deme 1]
	free(IDind_demes); // [0; 1; 2; 0; 1; 2]
	free(IDind_metapop); // [0; 1; 2; 3; 4; 5]
	free(maleAlloc_metapop); // [0.6; 0.7; 0.5; 0.5; 0.6; 0.9]
	free(femaleAlloc_metapop); // [0.4; 0.3; 0.5; 0.5; 0.4; 0.1]
}


void migration_step(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int extinctionStatus[], const int nImmigrants[]){
	// deals with migration : 
	//	1. emigrant demes are sampled according to their fecundity
	//	2. migrants are chosen uniformly within the sampled emigrant deme, and then cloned to the dest. deme
	int nEmigrantDemes = nDemes - 1; // number of emigrant demes = number of total demes, minus the receiving deme, minus the extincted demes
	
	int deme_tmp = 0;
	// remove the extincted demes
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if( extinctionStatus[deme_tmp]==1 ){
			nEmigrantDemes -= 1;
		}
	}
	
	if(nEmigrantDemes<=0){
		printf("There is no deme able to send migrants\n");
		exit(1);
	}
	
	// loop over receiving demes	
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
// CR		printf("\nCurrent deme: %d\tStatus: %d\tnIndividuals: %d\tnImmigrants: %d\n", deme_tmp, extinctionStatus[deme_tmp], newPopulation[deme_tmp].nIndividus, nImmigrants[deme_tmp]);
		if(extinctionStatus[deme_tmp]==0){ // only non extincted demes receive migrants
			// dealing with local seeds
			int nImmigrants_in_deme = 0; // number of coming seeds, i.e, number of local seeds to replace
			nImmigrants_in_deme = nImmigrants[deme_tmp];
			
			if(nImmigrants_in_deme>0){
				int* all_local_seeds = NULL; // list of seeds that can be replaced by migrants
				all_local_seeds = malloc(newPopulation[deme_tmp].nIndividus*sizeof(int));
				
				int* replaced_seeds = NULL; // list of local seeds that are replaced by migrants
				replaced_seeds = malloc(nImmigrants_in_deme*sizeof(int));
				
				if(all_local_seeds == NULL || replaced_seeds == NULL){
					exit(0);
				}
				
				int i = 0;
				for(i=0; i<newPopulation[deme_tmp].nIndividus; i++){
					all_local_seeds[i] = i;
				}
				
				gsl_ran_choose(r, replaced_seeds, nImmigrants_in_deme, all_local_seeds, newPopulation[deme_tmp].nIndividus, sizeof(int)); // list of locally replaced individuals in replaced_seeds (of length nImmigrants_in_deme)
/*				int pouet=0;
				printf("Individuals that can be replaced: ");
				for(pouet=0; pouet<newPopulation[deme_tmp].nIndividus; pouet++){
					printf("%d\t", all_local_seeds[pouet]);
				}
				printf("\n");
				
				printf("Replaced individuals: ");
				for(pouet=0; pouet<nImmigrants_in_deme; pouet++){
					printf("%d\t", replaced_seeds[pouet]);
				}
				printf("\n");
*/
			
				// dealing with emigrant demes
				double* emigrantDemes = NULL;
				double* fertilityPerDemes = NULL;
				emigrantDemes = malloc(nEmigrantDemes*sizeof(double)); // ID of demes allowing to send migrants
				fertilityPerDemes = malloc(nEmigrantDemes*sizeof(double)); // fertility of demes allowing to send migrants
				
				if( emigrantDemes == NULL || fertilityPerDemes == NULL ){
					exit(0);
				}
				
				int cnt = 0;
				for(i=0; i<nDemes; i++){
					if( extinctionStatus[i] == 0 ){
						if( i!=deme_tmp ){
							emigrantDemes[cnt] = (double)i;
							fertilityPerDemes[cnt]= population[i].fertility;
							cnt += 1;
						}
					}
				}
				
				double* sampled_emigrant_demes = NULL;
				sampled_emigrant_demes = malloc(nImmigrants_in_deme*sizeof(double));
				
				weightedSample(r, emigrantDemes, fertilityPerDemes, sampled_emigrant_demes, nEmigrantDemes, nImmigrants_in_deme); // sampled_emigrant_demes (of length nImmigrants_in_deme) contains the list of demes sending migrants
				
				// loop over migrants
				int migrant_tmp = 0;
				for(migrant_tmp=0; migrant_tmp<nImmigrants_in_deme; migrant_tmp++){
					// copy/paste version of the migration algorithm
					int copied_seed_tmp = 0; // position of the copied individual within its original deme ( sampled_emigrant_demes[migrant_tmp])
					copied_seed_tmp = gsl_ran_flat(r, 0, newPopulation[ (int)sampled_emigrant_demes[migrant_tmp] ].nIndividus); // randomly choose the emigrant individual from the emigrant deme
					
					int replaced_seed_tmp = 0; // position of the replaced individual within the current deme (deme_tmp)
					replaced_seed_tmp = replaced_seeds[migrant_tmp];
					
					// neutral alleles
					int ntrlLoci_tmp = 0;
					for(ntrlLoci_tmp=0; ntrlLoci_tmp<nNtrlLoci; ntrlLoci_tmp++){
						int substitute_pos = 0; // position of the allele replacing the local allele
						substitute_pos = (2*nNtrlLoci*copied_seed_tmp) + (2*ntrlLoci_tmp);
						
						int replaced_pos = 0; // position of the local allele
						replaced_pos = (2*nNtrlLoci*replaced_seed_tmp) + (2*ntrlLoci_tmp);
						
						int allele = 0;
						for(allele=0; allele<2; allele++){
							newPopulation[deme_tmp].ntrlLoci[replaced_pos + allele] = newPopulation[ (int)sampled_emigrant_demes[migrant_tmp] ].ntrlLoci[substitute_pos + allele]; // allele 0
						}
					}
					
					// quantitative loci
					newPopulation[deme_tmp].femaleAllocation[replaced_seed_tmp] = 0.0;
					int quantiLoci_tmp = 0;
					for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){
						int substitute_pos = 0; // position of the allele replacing the local allele
						substitute_pos = (2*nQuantiLoci*copied_seed_tmp) + (2*quantiLoci_tmp);
						
						int replaced_pos = 0; // position of the local allele
						replaced_pos = (2*nQuantiLoci*replaced_seed_tmp) + (2*quantiLoci_tmp);
						
						int allele = 0;
						for(allele=0; allele<2; allele++){
							newPopulation[deme_tmp].quantiLoci[replaced_pos + allele] = newPopulation[ (int)sampled_emigrant_demes[migrant_tmp] ].quantiLoci[substitute_pos + allele]; // allele 0
							newPopulation[deme_tmp].femaleAllocation[replaced_seed_tmp] += newPopulation[deme_tmp].quantiLoci[replaced_pos + allele];
						}
					}
					newPopulation[deme_tmp].maleAllocation[replaced_seed_tmp] = 1 - newPopulation[deme_tmp].femaleAllocation[replaced_seed_tmp];
				}
				free(replaced_seeds);
				free(emigrantDemes);
				free(fertilityPerDemes);
			} // end of the condition on the number of migrants>0
		} // end of the condition on the extinction status
	} // end of loop over receiving demes
}

void migration_step_v2(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int extinctionStatus[], const int nImmigrants[], const double selfingRate, const double dispersal){
	// deals with migration : 
	//	1. the migrants are produced not by cloning, but by sampling of mothers (z) and the fathers ((1-z)/(nDemes-1))
	int deme_tmp = 0;
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==0){
			int nMigrants = nImmigrants[deme_tmp];
			if(nMigrants>0){
				int local_deme = deme_tmp;
				
				double* deme_target_mother = NULL; // deme's IDs of the k mothers
				double* deme_target_father = NULL; // deme's IDs of the k fathers
				int* ind_target_mother = NULL; // ind's IDs of the k mothers
				int* ind_target_father = NULL; // ind's IDs of the k fathers

				deme_target_mother = malloc(nMigrants*sizeof(double)); // [deme #123; deme #412]
				deme_target_father = malloc(nMigrants*sizeof(double)); // [deme #3; deme #671]
				ind_target_mother = malloc(nMigrants*sizeof(int)); // [individual #20; individual #12]
				ind_target_father = malloc(nMigrants*sizeof(int)); // [individual #0; individual #2]
				
				if(deme_target_mother==NULL || deme_target_father==NULL || ind_target_mother==NULL || ind_target_father==NULL){
					exit(0);
				}
				
				// entries:
				//	r
				//	number of parents to sample
				//	the deme's IDs where the sampled father come from
				//	the deme's IDs where the sampled mother come from
				//	the individual's IDs within the sampled deme of the sampled father 
				//	the individual's IDs within the sampled deme of the sampled mother 
				//	the parental metapopulation
				//	number of demes in the metapopulation
				//	extinction status
				//	selfing rate
				//	pollen dispersal rate
				//	local deme where the dispersed seed will arrived
				//	model of dispersal; 0 = migrant pool model 1 = propagule pool model
				get_parent_for_dispersal(r, nMigrants, deme_target_mother, deme_target_father, ind_target_mother, ind_target_father, population, nDemes, extinctionStatus, selfingRate, dispersal, local_deme, 0);

				// making the seeds
				int ind_tmp = 0;
				for(ind_tmp=0; ind_tmp<nMigrants; ind_tmp++){
					int ind_to_replace = 0;
					ind_to_replace = gsl_rng_uniform_int(r, newPopulation[deme_tmp].nIndividus); // replace one random individual
					// neutral loci
					int ntrlLoci_tmp = 0;
					for(ntrlLoci_tmp=0; ntrlLoci_tmp<nNtrlLoci; ntrlLoci_tmp++){ // loop along the neutral loci
						// mother
						int pos_seed = 0; // position within the receiving deme
						int pos_parent = 0; // position within the donor deme
						pos_seed = 2*nNtrlLoci*ind_to_replace + 2*ntrlLoci_tmp + 0;
						pos_parent = 2*nNtrlLoci*ind_target_mother[ind_tmp] + 2*ntrlLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
						newPopulation[deme_tmp].ntrlLoci[pos_seed] = population[ (int)deme_target_mother[ind_tmp] ].ntrlLoci[pos_parent];
						
						// father
						pos_seed += 1;
						pos_parent = 2*nNtrlLoci*ind_target_father[ind_tmp] + 2*ntrlLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
						
						newPopulation[deme_tmp].ntrlLoci[pos_seed] = population[ (int)deme_target_father[ind_tmp] ].ntrlLoci[pos_parent];
					} // end of loop for neutral loci
					
					// quantitative loci
					int quantiLoci_tmp = 0;
					newPopulation[deme_tmp].femaleAllocation[ind_to_replace] = 0.0;
					for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){ // loop along the quantitative loci
						// mother
						int pos_seed = 0; // position within the receiving deme
						int pos_parent = 0; // position within the donor deme
						pos_seed = 2*nQuantiLoci*ind_to_replace + 2*quantiLoci_tmp + 0;
						pos_parent = 2*nQuantiLoci*ind_target_mother[ind_tmp] + 2*quantiLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
						
						newPopulation[deme_tmp].quantiLoci[pos_seed] = population[ (int)deme_target_mother[ind_tmp] ].quantiLoci[pos_parent];
						
						newPopulation[deme_tmp].femaleAllocation[ind_to_replace] += newPopulation[deme_tmp].quantiLoci[pos_seed];
						
						// father
						pos_seed += 1;
						pos_parent = 2*nQuantiLoci*ind_target_father[ind_tmp] + 2*quantiLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
						
						newPopulation[deme_tmp].quantiLoci[pos_seed] = population[ (int)deme_target_father[ind_tmp] ].quantiLoci[pos_parent];
						newPopulation[deme_tmp].femaleAllocation[ind_to_replace] += newPopulation[deme_tmp].quantiLoci[pos_seed];
						
					} // end of loop for quantitative loci
					// male allocation
					newPopulation[deme_tmp].maleAllocation[ind_to_replace] = 1-newPopulation[deme_tmp].femaleAllocation[ind_to_replace];
					
				} // end of loop to make seeds
			free(deme_target_mother);
			free(deme_target_father);
			free(ind_target_mother);
			free(ind_target_father);
			} // end of condition on the number of migrants
		} // end of condition on extinction
	} // end of loop over demes
}

void recolonization_step(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int extinctionStatus[], const int recolonization, const int colonizationModel, const double selfingRate, const int maxIndPerDem){
	// Deals with recolonization
	//	1. colonisers are produced by cloning random guys in the metapopulation according to the migrant-pool / propagule-pool model
	//	2. doesn't directly deals with pollen disperal
	double* list_demes = NULL; // contains all demes
	double* weight_demes = NULL; // weight of 0 for extincted, deme's fecundity for non extincted
	list_demes = malloc(nDemes*sizeof(double));
	weight_demes = malloc(nDemes*sizeof(double));
	
	int deme_tmp = 0;
	
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		list_demes[deme_tmp] = deme_tmp;
		if(extinctionStatus[deme_tmp]==0){
			weight_demes[deme_tmp] = population[deme_tmp].fertility;
		}else{
			weight_demes[deme_tmp] = 0;
		}
	}
	
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==1){
			double* deme_target_donor = NULL; // deme's IDs of the k donors 
			deme_target_donor = malloc(recolonization*sizeof(double)); // [deme #123; deme #412]
			
			if(deme_target_donor==NULL){
				exit(0);
			}
			
			// void weightedSample(gsl_rng* r, const double* liste, const double* weights, double* target, const int sizeOfListe, const int nTrials){
			weightedSample(r, list_demes, weight_demes, deme_target_donor, nDemes, recolonization);
			if(colonizationModel==1){ // if propagule pool model
				int i=0;
				for(i=0; i<recolonization; i++){
					deme_target_donor[i] = deme_target_donor[0];
				}
			}
			
			// making the seeds
			int ind_tmp = 0;
			for(ind_tmp=0; ind_tmp<recolonization; ind_tmp++){
				int copied_seed_tmp = 0; // position of the copied individual within its original deme ( sampled_emigrant_demes[migrant_tmp])
				copied_seed_tmp = gsl_ran_flat(r, 0, newPopulation[ (int)deme_target_donor[ind_tmp] ].nIndividus); // randomly choose the emigrant individual from the emigrant deme
				
				// neutral loci
				int ntrlLoci_tmp = 0;
				for(ntrlLoci_tmp=0; ntrlLoci_tmp<nNtrlLoci; ntrlLoci_tmp++){ // loop along the neutral loci
					// mother
					int pos_seed = 0; // position within the receiving deme
					int pos_parent = 0; // position within the donor deme
					pos_seed = 2*nNtrlLoci*ind_tmp + 2*ntrlLoci_tmp + 0;
					pos_parent = 2*nNtrlLoci*copied_seed_tmp + 2*ntrlLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					newPopulation[deme_tmp].ntrlLoci[pos_seed] = population[ (int)deme_target_donor[ind_tmp] ].ntrlLoci[pos_parent];
					
					// father
					pos_seed += 1;
					pos_parent = 2*nNtrlLoci*copied_seed_tmp + 2*ntrlLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					newPopulation[deme_tmp].ntrlLoci[pos_seed] = population[ (int)deme_target_donor[ind_tmp] ].ntrlLoci[pos_parent];
				} // end of loop for neutral loci
				
				// quantitative loci
				int quantiLoci_tmp = 0;
				newPopulation[deme_tmp].femaleAllocation[ind_tmp] = 0.0;
				for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){ // loop along the quantitative loci
					// mother
					int pos_seed = 0; // position within the receiving deme
					int pos_parent = 0; // position within the donor deme
					pos_seed = 2*nQuantiLoci*ind_tmp + 2*quantiLoci_tmp + 0;
					pos_parent = 2*nQuantiLoci*copied_seed_tmp + 2*quantiLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].quantiLoci[pos_seed] = population[ (int)deme_target_donor[ind_tmp] ].quantiLoci[pos_parent];
					
					newPopulation[deme_tmp].femaleAllocation[ind_tmp] += newPopulation[deme_tmp].quantiLoci[pos_seed];
					
					// father
					pos_seed += 1;
					pos_parent = 2*nQuantiLoci*copied_seed_tmp + 2*quantiLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].quantiLoci[pos_seed] = population[ (int)deme_target_donor[ind_tmp] ].quantiLoci[pos_parent];
					newPopulation[deme_tmp].femaleAllocation[ind_tmp] += newPopulation[deme_tmp].quantiLoci[pos_seed];
					
				} // end of loop for quantitative loci
				// male allocation
				newPopulation[deme_tmp].maleAllocation[ind_tmp] = 1-newPopulation[deme_tmp].femaleAllocation[ind_tmp];
				
			} // end of loop to make seeds
			
			// second round of reproduction to fill the deme
			// // get the second mothers
			const int seeds_second = maxIndPerDem;
			double* mothers_second = NULL;
			mothers_second = malloc(seeds_second * sizeof(double));
			
			double* female_fitness_second = NULL;
			female_fitness_second = malloc(recolonization * sizeof(double));
			
			double* females_second = NULL;
			females_second= malloc(recolonization * sizeof(double));
			
			for(ind_tmp=0; ind_tmp<recolonization; ind_tmp++){
				females_second[ind_tmp] = ind_tmp;
				female_fitness_second[ind_tmp] = newPopulation[deme_tmp].femaleAllocation[ind_tmp];
			}
			
			weightedSample(r, females_second, female_fitness_second, mothers_second, recolonization, seeds_second);
			
			// // get the second fathers
			double* fathers_second = NULL;
			fathers_second = malloc(seeds_second * sizeof(double));
			
			double* male_fitness_second = NULL;
			male_fitness_second = malloc(recolonization * sizeof(double));
			
			double* males_second = NULL;
			males_second= malloc(recolonization * sizeof(double));
			
			for(ind_tmp=0; ind_tmp<recolonization; ind_tmp++){
				males_second[ind_tmp] = ind_tmp;
				male_fitness_second[ind_tmp] = newPopulation[deme_tmp].maleAllocation[ind_tmp];
			}
			
			weightedSample(r, males_second, male_fitness_second, fathers_second, recolonization, seeds_second);
			
			// // selfing 
			for(ind_tmp=0; ind_tmp<seeds_second; ind_tmp++){
				int autofec = 0;
				autofec = gsl_ran_binomial(r, selfingRate, 1);
				if(autofec==1){
					fathers_second[ind_tmp] = mothers_second[ind_tmp];
				}
			}
			
			// // produce babies
			//int nIndividus;
			//long* ntrlLoci;	// contains the neutral alleles for nIndividus (number of individuals within the deme) x 2 (diploids) x nNtlrLoci (number of neutral loci)
			//double* quantiLoci; // contains the allelic effects for nIndividus (number of individuals within the deme) x 2 (diploids) x nQuantiLoci (number of quantitative loci)
			//double* femaleAllocation; // =sum of the allelic effects over the quantitative loci
			//double* maleAllocation; // =(1 - femaleAllocation)
			//double fertility; // sum of individual female allocations within the deme
			// 1) neutral loci
			long* ntrlLoci_second = NULL;
			ntrlLoci_second = malloc(2*nNtrlLoci*seeds_second * sizeof(long));
			
			for(ind_tmp=0; ind_tmp<seeds_second; ind_tmp++){
				int ntrlLoci_tmp = 0;
				for(ntrlLoci_tmp=0; ntrlLoci_tmp<nNtrlLoci; ntrlLoci_tmp++){
					// from mother
					int pos_seed = 0;
					pos_seed = (2*nNtrlLoci*ind_tmp) + (2*ntrlLoci_tmp) + 0;
					int pos_parent = 0;
					pos_parent = (2*nNtrlLoci*mothers_second[ind_tmp]) + (2*ntrlLoci_tmp) + gsl_ran_binomial(r, 0.5, 1);
					
					ntrlLoci_second[pos_seed] = newPopulation[ deme_tmp ].ntrlLoci[pos_parent];
					
					// from father 
					pos_seed += 1;
					pos_parent = (2*nNtrlLoci*fathers_second[ind_tmp]) + (2*ntrlLoci_tmp) + gsl_ran_binomial(r, 0.5, 1);
					
					ntrlLoci_second[pos_seed] = newPopulation[ deme_tmp ].ntrlLoci[pos_parent];
				}
			}
			
			// 2) quantitative loci
			double* quantiLoci_second = NULL;
			quantiLoci_second = malloc(2*nQuantiLoci*seeds_second * sizeof(double));
			
			for(ind_tmp=0; ind_tmp<seeds_second; ind_tmp++){
				int quantiLoci_tmp = 0;
				for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){
					// from mother
					int pos_seed = 0;
					pos_seed = (2*nQuantiLoci*ind_tmp) + (2*quantiLoci_tmp) + 0;
					int pos_parent = 0;
					pos_parent = (2*nQuantiLoci*mothers_second[ind_tmp]) + (2*quantiLoci_tmp) + gsl_ran_binomial(r, 0.5, 1);
					
					quantiLoci_second[pos_seed] = newPopulation[ deme_tmp ].quantiLoci[pos_parent];
					
					// from father 
					pos_seed += 1;
					pos_parent = (2*nQuantiLoci*fathers_second[ind_tmp]) + (2*quantiLoci_tmp) + gsl_ran_binomial(r, 0.5, 1);
					
					quantiLoci_second[pos_seed] = newPopulation[ deme_tmp ].quantiLoci[pos_parent];
				}
			}
			
			// 3) sex allocation
			double* femaleAllocation_second = NULL;
			femaleAllocation_second = malloc(seeds_second * sizeof(double));
			double* maleAllocation_second = NULL;
			maleAllocation_second = malloc(seeds_second * sizeof(double));
			double fertility_second = 0.0;
			
			ind_tmp=0;
			for(ind_tmp=0; ind_tmp<seeds_second; ind_tmp++){
				femaleAllocation_second[ind_tmp] = 0.0;
				int quantiLoci_tmp = 0;
				for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){
					femaleAllocation_second[ind_tmp] += quantiLoci_second[ (2*nQuantiLoci*ind_tmp) + (2*quantiLoci_tmp) + 0 ];
					femaleAllocation_second[ind_tmp] += quantiLoci_second[ (2*nQuantiLoci*ind_tmp) + (2*quantiLoci_tmp) + 1 ];

					fertility_second += femaleAllocation_second[ind_tmp];
				}
				maleAllocation_second[ind_tmp] = 1-femaleAllocation_second[ind_tmp];
			}
			
			// Replace the individuals
			newPopulation[deme_tmp].nIndividus = seeds_second;
			newPopulation[deme_tmp].fertility = fertility_second;

			free(newPopulation[deme_tmp].femaleAllocation);			
			newPopulation[deme_tmp].femaleAllocation = malloc(seeds_second * sizeof(double));
			int i=0;
			for(i=0; i<seeds_second; i++){
				newPopulation[deme_tmp].femaleAllocation[i] = femaleAllocation_second[i];
			}
			
			free(newPopulation[deme_tmp].maleAllocation);			
			newPopulation[deme_tmp].maleAllocation = malloc(seeds_second * sizeof(double));
			for(i=0; i<seeds_second; i++){
				newPopulation[deme_tmp].maleAllocation[i] = maleAllocation_second[i];
			}
			
			free(newPopulation[deme_tmp].ntrlLoci);			
			newPopulation[deme_tmp].ntrlLoci = malloc(2*nNtrlLoci*seeds_second * sizeof(long));
			for(i=0; i<(2*seeds_second*nNtrlLoci); i++){
				newPopulation[deme_tmp].ntrlLoci[i] = ntrlLoci_second[i];
			}
			
			free(newPopulation[deme_tmp].quantiLoci);			
			newPopulation[deme_tmp].quantiLoci = malloc(2*nQuantiLoci*seeds_second * sizeof(long));
			for(i=0; i<(2*seeds_second*nQuantiLoci); i++){
				newPopulation[deme_tmp].quantiLoci[i] = quantiLoci_second[i];
			}
			
			
			free(mothers_second);
			free(female_fitness_second);
			free(females_second);
			free(fathers_second);
			free(male_fitness_second);
			free(males_second);
			free(ntrlLoci_second);
			free(quantiLoci_second);
			free(femaleAllocation_second);
			free(maleAllocation_second);
			
			free(deme_target_donor);
		} // end of condition on the extinction
	} // end of loop over demes
	free(list_demes);
	free(weight_demes);
} // end of function


void recolonization_step_v2(gsl_rng* r, const Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int extinctionStatus[], const int recolonization, const int colonizationModel, const double dispersal, const double selfingRate, const int maxIndPerDem){
	// Deals with recolonization
	//	1. colonisers are produced by sampling mothers and fathers from the whole metapopulation, according to (1-z) and (z). Not by cloning.
	int deme_tmp = 0;
	
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==1){
			int local_deme = deme_tmp;
			
			double* deme_target_mother = NULL; // deme's IDs of the k mothers
			double* deme_target_father = NULL; // deme's IDs of the k fathers
			int* ind_target_mother = NULL; // ind's IDs of the k mothers
			int* ind_target_father = NULL; // ind's IDs of the k fathers

			deme_target_mother = malloc(recolonization*sizeof(double)); // [deme #123; deme #412]
			deme_target_father = malloc(recolonization*sizeof(double)); // [deme #3; deme #671]
			ind_target_mother = malloc(recolonization*sizeof(int)); // [individual #20; individual #12]
			ind_target_father = malloc(recolonization*sizeof(int)); // [individual #0; individual #2]
			
			if(deme_target_mother==NULL || deme_target_father==NULL || ind_target_mother==NULL || ind_target_father==NULL){
				exit(0);
			}
			
			// entries:
			//	r
			//	number of parents to sample
			//	the deme's IDs where the sampled father come from
			//	the deme's IDs where the sampled mother come from
			//	the individual's IDs within the sampled deme of the sampled father 
			//	the individual's IDs within the sampled deme of the sampled mother 
			//	the parental metapopulation
			//	number of demes in the metapopulation
			//	extinction status
			//	selfing rate
			//	pollen dispersal rate
			//	local deme where the dispersed seed will arrived
			//	model of dispersal; 0 = migrant pool model 1 = propagule pool model
			get_parent_for_dispersal(r, recolonization, deme_target_mother, deme_target_father, ind_target_mother, ind_target_father, population, nDemes, extinctionStatus, selfingRate, dispersal, local_deme, colonizationModel);

			// making the seeds
			int ind_tmp = 0;
			for(ind_tmp=0; ind_tmp<recolonization; ind_tmp++){
				// neutral loci
				int ntrlLoci_tmp = 0;
				for(ntrlLoci_tmp=0; ntrlLoci_tmp<nNtrlLoci; ntrlLoci_tmp++){ // loop along the neutral loci
					// mother
					int pos_seed = 0; // position within the receiving deme
					int pos_parent = 0; // position within the donor deme
					pos_seed = 2*nNtrlLoci*ind_tmp + 2*ntrlLoci_tmp + 0;
					pos_parent = 2*nNtrlLoci*ind_target_mother[ind_tmp] + 2*ntrlLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					newPopulation[deme_tmp].ntrlLoci[pos_seed] = population[ (int)deme_target_mother[ind_tmp] ].ntrlLoci[pos_parent];
					
					// father
					pos_seed += 1;
					pos_parent = 2*nNtrlLoci*ind_target_father[ind_tmp] + 2*ntrlLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].ntrlLoci[pos_seed] = population[ (int)deme_target_father[ind_tmp] ].ntrlLoci[pos_parent];
				} // end of loop for neutral loci
				
				// quantitative loci
				int quantiLoci_tmp = 0;
				newPopulation[deme_tmp].femaleAllocation[ind_tmp] = 0.0;
				for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){ // loop along the quantitative loci
					// mother
					int pos_seed = 0; // position within the receiving deme
					int pos_parent = 0; // position within the donor deme
					pos_seed = 2*nQuantiLoci*ind_tmp + 2*quantiLoci_tmp + 0;
					pos_parent = 2*nQuantiLoci*ind_target_mother[ind_tmp] + 2*quantiLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].quantiLoci[pos_seed] = population[ (int)deme_target_mother[ind_tmp] ].quantiLoci[pos_parent];
					
					newPopulation[deme_tmp].femaleAllocation[ind_tmp] += newPopulation[deme_tmp].quantiLoci[pos_seed];
					
					// father
					pos_seed += 1;
					pos_parent = 2*nQuantiLoci*ind_target_father[ind_tmp] + 2*quantiLoci_tmp + gsl_ran_binomial(r, 0.5, 1);
					
					newPopulation[deme_tmp].quantiLoci[pos_seed] = population[ (int)deme_target_father[ind_tmp] ].quantiLoci[pos_parent];
					newPopulation[deme_tmp].femaleAllocation[ind_tmp] += newPopulation[deme_tmp].quantiLoci[pos_seed];
					
				} // end of loop for quantitative loci
				// male allocation
				newPopulation[deme_tmp].maleAllocation[ind_tmp] = 1-newPopulation[deme_tmp].femaleAllocation[ind_tmp];
				
			} // end of loop to make seeds
			
			// second round of reproduction to fill the deme
			// // get the second mothers
			const int seeds_second = maxIndPerDem;
			double* mothers_second = NULL;
			mothers_second = malloc(seeds_second * sizeof(double));
			
			double* female_fitness_second = NULL;
			female_fitness_second = malloc(recolonization * sizeof(double));
			
			double* females_second = NULL;
			females_second= malloc(recolonization * sizeof(double));
			
			for(ind_tmp=0; ind_tmp<recolonization; ind_tmp++){
				females_second[ind_tmp] = ind_tmp;
				female_fitness_second[ind_tmp] = newPopulation[deme_tmp].femaleAllocation[ind_tmp];
			}
			
			weightedSample(r, females_second, female_fitness_second, mothers_second, recolonization, seeds_second);
			
			// // get the second fathers
			double* fathers_second = NULL;
			fathers_second = malloc(seeds_second * sizeof(double));
			
			double* male_fitness_second = NULL;
			male_fitness_second = malloc(recolonization * sizeof(double));
			
			double* males_second = NULL;
			males_second= malloc(recolonization * sizeof(double));
			
			for(ind_tmp=0; ind_tmp<recolonization; ind_tmp++){
				males_second[ind_tmp] = ind_tmp;
				male_fitness_second[ind_tmp] = newPopulation[deme_tmp].maleAllocation[ind_tmp];
			}
			
			weightedSample(r, males_second, male_fitness_second, fathers_second, recolonization, seeds_second);
			
			// // selfing 
			for(ind_tmp=0; ind_tmp<seeds_second; ind_tmp++){
				int autofec = 0;
				autofec = gsl_ran_binomial(r, selfingRate, 1);
				if(autofec==1){
					fathers_second[ind_tmp] = mothers_second[ind_tmp];
				}
			}
			
			// // produce babies
			//int nIndividus;
			//long* ntrlLoci;	// contains the neutral alleles for nIndividus (number of individuals within the deme) x 2 (diploids) x nNtlrLoci (number of neutral loci)
			//double* quantiLoci; // contains the allelic effects for nIndividus (number of individuals within the deme) x 2 (diploids) x nQuantiLoci (number of quantitative loci)
			//double* femaleAllocation; // =sum of the allelic effects over the quantitative loci
			//double* maleAllocation; // =(1 - femaleAllocation)
			//double fertility; // sum of individual female allocations within the deme
			// 1) neutral loci
			long* ntrlLoci_second = NULL;
			ntrlLoci_second = malloc(2*nNtrlLoci*seeds_second * sizeof(long));
			
			for(ind_tmp=0; ind_tmp<seeds_second; ind_tmp++){
				int ntrlLoci_tmp = 0;
				for(ntrlLoci_tmp=0; ntrlLoci_tmp<nNtrlLoci; ntrlLoci_tmp++){
					// from mother
					int pos_seed = 0;
					pos_seed = (2*nNtrlLoci*ind_tmp) + (2*ntrlLoci_tmp) + 0;
					int pos_parent = 0;
					pos_parent = (2*nNtrlLoci*mothers_second[ind_tmp]) + (2*ntrlLoci_tmp) + gsl_ran_binomial(r, 0.5, 1);
					
					ntrlLoci_second[pos_seed] = newPopulation[ deme_tmp ].ntrlLoci[pos_parent];
					
					// from father 
					pos_seed += 1;
					pos_parent = (2*nNtrlLoci*fathers_second[ind_tmp]) + (2*ntrlLoci_tmp) + gsl_ran_binomial(r, 0.5, 1);
					
					ntrlLoci_second[pos_seed] = newPopulation[ deme_tmp ].ntrlLoci[pos_parent];
				}
			}
			
			// 2) quantitative loci
			double* quantiLoci_second = NULL;
			quantiLoci_second = malloc(2*nQuantiLoci*seeds_second * sizeof(double));
			
			for(ind_tmp=0; ind_tmp<seeds_second; ind_tmp++){
				int quantiLoci_tmp = 0;
				for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){
					// from mother
					int pos_seed = 0;
					pos_seed = (2*nQuantiLoci*ind_tmp) + (2*quantiLoci_tmp) + 0;
					int pos_parent = 0;
					pos_parent = (2*nQuantiLoci*mothers_second[ind_tmp]) + (2*quantiLoci_tmp) + gsl_ran_binomial(r, 0.5, 1);
					
					quantiLoci_second[pos_seed] = newPopulation[ deme_tmp ].quantiLoci[pos_parent];
					
					// from father 
					pos_seed += 1;
					pos_parent = (2*nQuantiLoci*fathers_second[ind_tmp]) + (2*quantiLoci_tmp) + gsl_ran_binomial(r, 0.5, 1);
					
					quantiLoci_second[pos_seed] = newPopulation[ deme_tmp ].quantiLoci[pos_parent];
				}
			}
			
			// 3) sex allocation
			double* femaleAllocation_second = NULL;
			femaleAllocation_second = malloc(seeds_second * sizeof(double));
			double* maleAllocation_second = NULL;
			maleAllocation_second = malloc(seeds_second * sizeof(double));
			double fertility_second = 0.0;
			
			ind_tmp=0;
			for(ind_tmp=0; ind_tmp<seeds_second; ind_tmp++){
				femaleAllocation_second[ind_tmp] = 0.0;
				int quantiLoci_tmp = 0;
				for(quantiLoci_tmp=0; quantiLoci_tmp<nQuantiLoci; quantiLoci_tmp++){
					femaleAllocation_second[ind_tmp] += quantiLoci_second[ (2*nQuantiLoci*ind_tmp) + (2*quantiLoci_tmp) + 0 ];
					femaleAllocation_second[ind_tmp] += quantiLoci_second[ (2*nQuantiLoci*ind_tmp) + (2*quantiLoci_tmp) + 1 ];

					fertility_second += femaleAllocation_second[ind_tmp];
				}
				maleAllocation_second[ind_tmp] = 1-femaleAllocation_second[ind_tmp];
			}
			
			// Replace the individuals
			newPopulation[deme_tmp].nIndividus = seeds_second;
			newPopulation[deme_tmp].fertility = fertility_second;

			free(newPopulation[deme_tmp].femaleAllocation);			
			newPopulation[deme_tmp].femaleAllocation = malloc(seeds_second * sizeof(double));
			int i=0;
			for(i=0; i<seeds_second; i++){
				newPopulation[deme_tmp].femaleAllocation[i] = femaleAllocation_second[i];
			}
			
			free(newPopulation[deme_tmp].maleAllocation);			
			newPopulation[deme_tmp].maleAllocation = malloc(seeds_second * sizeof(double));
			for(i=0; i<seeds_second; i++){
				newPopulation[deme_tmp].maleAllocation[i] = maleAllocation_second[i];
			}
			
			free(newPopulation[deme_tmp].ntrlLoci);			
			newPopulation[deme_tmp].ntrlLoci = malloc(2*nNtrlLoci*seeds_second * sizeof(long));
			for(i=0; i<(2*seeds_second*nNtrlLoci); i++){
				newPopulation[deme_tmp].ntrlLoci[i] = ntrlLoci_second[i];
			}
			
			free(newPopulation[deme_tmp].quantiLoci);			
			newPopulation[deme_tmp].quantiLoci = malloc(2*nQuantiLoci*seeds_second * sizeof(long));
			for(i=0; i<(2*seeds_second*nQuantiLoci); i++){
				newPopulation[deme_tmp].quantiLoci[i] = quantiLoci_second[i];
			}
			
			
			free(mothers_second);
			free(female_fitness_second);
			free(females_second);
			free(fathers_second);
			free(male_fitness_second);
			free(males_second);
			free(ntrlLoci_second);
			free(quantiLoci_second);
			free(femaleAllocation_second);
			free(maleAllocation_second);
			
			free(deme_target_mother);
			free(deme_target_father);
			free(ind_target_mother);
			free(ind_target_father);
		} // end of condition on the extinction
	} // end of loop over demes
} // end of function

void weightedSample(gsl_rng* r, const double* liste, const double* weights, double* target, const int sizeOfListe, const int nTrials){
	// function that fills the vector 'target' of size 'nTrials' containing the weighted-sampled 'sizeOfListe' elements of the 'liste':
	// weightedSample(gsl_rng* r, {2, 4, 6, 8, 10}, {1.2, 0.6, 0.3, 0.15, 0.05}, target, 5, 20)
	// target = {6, 6, 6, 2, 2, 2, 2, 2, 8, 2, 8, 6, 4, 2, 2, 4, 4, 4, 4, 2}
	// but can also be used for boolean sampling (pile ou face) using:
	// weightedSample(gsl_rng* r, {0, 1}, {1, 1}, target, 2, 1)
	

	int i = 0;
	unsigned int* n = NULL;
	int* sampledListe = NULL;
	n = malloc(sizeOfListe * sizeof(double));	// will contain the number of succes after K nTrials for each of the sizeOfListe elements of liste
	sampledListe = malloc(nTrials * sizeof(int));	// if n={0, 3, 1, 1}, sampledListe={1, 1, 1, 4, 5}

	gsl_ran_multinomial(r, sizeOfListe, nTrials, weights, n);	// return in 'n' the number of success for the sizeOfListe elements of liste

	int nValues = 0;
	int tmp = 0;
	int tmp2 = 0;
	for(i=0; i<sizeOfListe; i++){ // loop along the list called 'n' resulting from gsl_ran_multinomial
		nValues = n[i];
		if(nValues != 0){
			tmp2 = 0;
			do{
				sampledListe[tmp] = i;
				tmp++;
				tmp2++;
			}while(tmp2 < nValues);
		}
	}

	// shuffle values of the sampledListe
	gsl_permutation* p = gsl_permutation_alloc (nTrials);
	gsl_permutation_init (p);
	gsl_ran_shuffle(r, p -> data, nTrials, sizeof(size_t));

	tmp = 0;
	for(i=0; i<nTrials; i++){
		tmp=gsl_permutation_get(p, i);
		target[i]=liste[sampledListe[tmp]];
	}
	gsl_permutation_free(p);
	free(n);
	free(sampledListe);
}

void libererMemoirePopulation(Deme* population, const int nDemes){
	// free memory taken by the population at the end of each generation.
	int i = 0;
	for(i=0; i<nDemes; i++){
		free(population[i].ntrlLoci);
		free(population[i].quantiLoci);
		free(population[i].femaleAllocation);
		free(population[i].maleAllocation);
	}
}

void writeNindividuals(const Deme* population, const int nDemes, const double extinction, const double migrationRate, const int seed){
	int i = 0;

	char nameOfFileNInd[100];
	char nameOfFileFemAlloc[100];
	sprintf(nameOfFileNInd, "indOverTime_e%f_i%f_seed%d.txt", extinction, migrationRate, seed);
	sprintf(nameOfFileFemAlloc, "femAllocOverTime_e%f_i%f_seed%d.txt", extinction, migrationRate, seed);

	FILE* fichierNInd = NULL;
	FILE* fichierFemAlloc = NULL;
	fichierNInd =  fopen(nameOfFileNInd, "a");
	fichierFemAlloc =  fopen(nameOfFileFemAlloc, "a");
	if(fichierNInd != NULL && fichierFemAlloc != NULL){
		for(i=0; i<nDemes; i++){
			fprintf(fichierNInd, "%d ", population[i].nIndividus);
			fprintf(fichierFemAlloc, "%f ", gsl_stats_mean(population[i].femaleAllocation, 1, population[i].nIndividus));
		}
		fprintf(fichierNInd, "\n");
		fclose(fichierNInd);
		fprintf(fichierFemAlloc, "\n");
		fclose(fichierFemAlloc);
	}
}

void statisticsPopulations(Deme* population, const int nDemes, const int maxIndPerDem, const int nQuantiLoci, const int fecundity, const double migrationRate, const double extinction, const double dispersal, const int recolonization, const int sexualSystem, const int seed, int time, const double selfingRate, const int colonizationModel, const double global_fst_cm, const double global_fst_coal, const double global_gpst, const double global_D, const double global_Fis, double avg_nc, double avg_Hs, double avg_Htot, double avg_Ho, const double global_Fst_anova, const double global_Fis_anova, const double global_Fit_anova, const double z_exp){
	// function that calculates the mean female allocation, its standard deviation and the percentage of cosexuals in the metapopulation
	// also computes FST_var = var(p)/(1-p) and FST_coal = (Htot - Hs) / Htot
	double fstValue = 0.0;
	double fstRoussetValue = 0.0;
	double fstValueDensity = 0.0;
	double meanAllocFemale = 0.0;
	double sdAllocFemale = 0.0;

	char nomFichierSortie[200];
	sprintf(nomFichierSortie, "output_%d.txt", seed);
	FILE* fichierSortie = NULL;

	fichierSortie = fopen(nomFichierSortie, "r");
	if(fichierSortie == NULL){
		fichierSortie = fopen(nomFichierSortie, "a");
		//fprintf(fichierSortie, "nDemes\tnIndMaxPerDeme\tNtot\tnQuantiLoci\tselfingRate\tfecundity\tmigRate\textRate\tcolonizationModel\trecolonization\tatGeneration\tsexSystem\tsexAvantage\tseed\tmeanFemAlloc\tsdFemAlloc\tmeanFemAllocCosexual\tsdFemAllocCosexual\tcosexualProportion\tobsFST_var\tobsFST_coal\tobsGST_p\tobsJostD\tobsFIS\texpFST_Nmax\texpFST_Nobs\texpFST_Rousset_Nmax\n");
		fprintf(fichierSortie, "nDemes\tnIndMaxPerDeme\tNtot\tnQuantiLoci\tselfingRate\tfecundity\tmigRate\textRate\tpollenDispersal\tcolonizationModel\trecolonization\tatGeneration\tsexSystem\tseed\tmeanFemAlloc\tsdFemAlloc\tobs_nc\tobs_Hs\tobs_Htot\tobs_Ho\tobsFST_var\tobsFST_coal\tobsGST_p\tobsJostD\tobsFIS\tobsFst_anova\tobsFis_anova\tobsFit_anova\texpFST_Nmax\texpFST_Nobs\texpFST_Rousset_Nmax\texp_femAlloc\n"); 
		fclose(fichierSortie);
	}else{
		fclose(fichierSortie);
	}
	
	fichierSortie = fopen(nomFichierSortie, "a");

	if(fichierSortie != NULL){	
		
		int nIndividusTotal = 0;
		int i = 0;
		for(i=0; i<nDemes; i++){
			nIndividusTotal += population[i].nIndividus;
		}	
	
		double* allocFemale = NULL; // female allocation in the whole metapopulation
		allocFemale = malloc(nIndividusTotal * sizeof(double));
		
		if(allocFemale==NULL){
			exit(0);
		}

		int cnt = 0;
		for(i=0; i<nDemes; i++){
			int j = 0;
			for(j=0; j<population[i].nIndividus; j++){
				allocFemale[cnt] = population[i].femaleAllocation[j];
				cnt += 1;
			}
		}

		meanAllocFemale = gsl_stats_mean(allocFemale, 1, nIndividusTotal);
		sdAllocFemale = gsl_stats_sd(allocFemale, 1, nIndividusTotal);
		
		fstValue = fstMullon(maxIndPerDem, extinction, recolonization, migrationRate); // expected fst assuming that all demes are full
		fstValueDensity = fstMullon((int) nIndividusTotal/(1.0*nDemes), extinction, recolonization, migrationRate); // expected fst, using the the average number of individuals in demes (for cases with high extinction, low fecundity)
		fstRoussetValue = fstRousset(maxIndPerDem, extinction, recolonization, migrationRate, colonizationModel);

		free(allocFemale);
		
		char colonizationModelTMP[50];
		
		if(colonizationModel == 0){
			sprintf(colonizationModelTMP, "migrationPool");
//			colonizationModelTMP = "migrationPool";
		}
		if(colonizationModel == 1){
			sprintf(colonizationModelTMP, "propagulePool");
//			colonizationModelTMP = "propagulePool";
		}

		fprintf(fichierSortie, "%d\t%d\t%d\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%s\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", nDemes, maxIndPerDem, nIndividusTotal, nQuantiLoci, selfingRate, fecundity, migrationRate, extinction, dispersal, colonizationModelTMP, recolonization, time, sexualSystem, seed, meanAllocFemale, sdAllocFemale, avg_nc, avg_Hs, avg_Htot, avg_Ho, global_fst_cm, global_fst_coal, global_gpst, global_D, global_Fis, global_Fst_anova, global_Fis_anova, global_Fit_anova, fstValue, fstValueDensity, fstRoussetValue, z_exp);
		fclose(fichierSortie);
	}
}

void checkCommandLine(int argc){
	if(argc != 18){
		printf("\n%sThe number of provided arguments is not correct.%s\nYou provided %s%d%s argument while %s19%s are expected:\n\t\
		%s1.%s  Number of demes (>0)\n\t\
		%s2.%s  Max number of individuals per deme (>0)\n\n\t\
		%s3.%s  Number of generations (>0)\n\n\t\
		%s4.%s  Number of neutral loci (>=0)\n\t\
		%s5.%s  Neutral mutation rate (in [0-1]))\n\t\
		%s6.%s  Number of quantitative loci (>0)\n\t\
		%s7.%s  Quantitative mutation rate (in [0-1])\n\n\t\
		%s8.%s  Max number of offsprings per hermaphrodite (>0)\n\n\t\
		%s9.%s Immigration rate (Poisson distributed; >=0)\n\t\
		%s10.%s Rate of pollen dispersal\n\n\t\
		%s11.%s Extinction rate (Binomialy distributed; in [0-1])\n\t\
		%s12.%s Number of individuals recolonizing an extincted deme (>0)\n\t\
		%s13.%s Colonization model, 'migration pool' (=0) or 'propagule pool' (=1) models\n\n\t\
		%s14.%s sexualSystem is equal to 0 if autosomal, equal to 1 if XY and equal to 2 if ZW\n\t\
		%s15.%s Selfing rate of hermaphrodites, fixed over time (in [0-1])\n\n\t\
		%s16.%s Frequency at which statistics are written in the outfile (>0)\n\n\t\
		%s17.%s Seed for the random generator (>0)\n\n\
		%s\tExample:%s ./quantiSex 100 100   100   10 0.0001 1 0.00001   100   1 0.05   0.1 1 1   0 0.1   20   123\n\n\
		version: %s\n\n\t\tdependencies: \t%s\n\n", KRED, STOP, KRED, argc-1, STOP, KRED, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KRED, STOP, VERSION, DEPENDENCY);

		exit(0);
	}
}

double fstMullon(const int maxIndPerDem, const double extinction, const int recolonization, const double migrationRate){
	double res = 0.0;
	double numQ = 0.0;
	double denomQ = 0.0;
//	double phi = 0.0;
	double qr = 0.0;
	double migrationProportion = migrationRate / maxIndPerDem;
	
	numQ = 1/(2.0 * maxIndPerDem) + extinction/(2.0 * recolonization) - extinction/(2.0 * recolonization * 2.0 * maxIndPerDem);
	denomQ = 1 - (1 - 1/(2.0 * maxIndPerDem)) * ((1 - migrationProportion)*(1 - migrationProportion) * (1 - extinction) + extinction * (1 - 1/(2.0 * recolonization)) * 1/(2.0 * recolonization -1));
//	phi = 1/(2.0 * recolonization -1);
	qr = numQ/denomQ;
	res = (qr - 1/(2.0 * maxIndPerDem)) * (2 * maxIndPerDem)/(2.0 * maxIndPerDem -1);

	return(res);
}

double fstRousset(const int maxIndPerDem, const double extinction, const int recolonization, const double migrationRate, const int colonizationModel){
	double res = 0.0;
	double numQ = 0.0;
	double denomQ = 0.0;
	double phi = 0.0;
	double qr = 0.0;
	double migrationProportion = migrationRate / maxIndPerDem;
	
	if(colonizationModel == 1){
		phi = pow((1 - migrationProportion), 2);
	}

	numQ = 1/(2.0 * maxIndPerDem) + extinction/(2.0 * recolonization) - extinction/(2.0 * recolonization * 2.0 * maxIndPerDem);
	denomQ = 1 - (1 - 1/(2.0 * maxIndPerDem)) * (pow((1 - migrationProportion), 2) * (1 - extinction) + extinction * phi * (1 - 1/(2.0 * recolonization)) * 1/(2.0 * recolonization -1));

	qr = numQ/denomQ;

	res = (qr - 1/(2.0 * maxIndPerDem)) * (2 * maxIndPerDem)/(2.0 * maxIndPerDem -1);

	return(res);
}

/*# diveRsityR
#!/usr/bin/env Rscript
./diveRsity.R input=nameOfGenePopFile output=nameOfROutputFile
library(diveRsity)
options(warn=-1)
for(i in commandArgs()){
	tmp = strsplit(i, "=")
	if(tmp[[1]][1] == "input"){input = tmp[[1]][2]}
	if(tmp[[1]][1] == "output"){output = tmp[[1]][2]}
}

a=diffCalc(input, fst=T, pairwise=F, outfile=output)

*/

void QoneQtwoQthree(const Deme* population, const int nDemes, const int nNtrlLoci, const int locus, double* target){
	int nTot = 0;
	int i = 0;
	int j = 0;
	long allele1 = 0;
	long allele2 = 0;
	double Q1 = 0.0; // proportion of homozygous individus: nHomozygous / nTot
	double Q2 = 0.0; // proportion of homozygous within demes. It's the average over demes, weighted by the number of pairs within demes
	double Q3 = 0.0; // proportion of homozygous among demes

	// contengency table of alleles for the whole metapop	
	double* cont_table_tot = NULL;
	cont_table_tot = malloc((MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1)*sizeof(double));
	
	// contengency table of alleles for each deme	
	double cont_table_demes[nDemes][MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1];
	
	if(cont_table_tot == NULL){
		exit(0);
	}

	for(i=0; i<nDemes; i++){
		for(j=0; j<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; j++){
			cont_table_demes[i][j] = 0.0;
		}
	}
	
	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		cont_table_tot[i] = 0;
	}	

	for(i=0; i<nDemes; i++){
		for(j=0; j<population[i].nIndividus; j++){
			nTot += 1;
			allele1 = population[i].ntrlLoci[2*nNtrlLoci*j + locus*2];
			allele2 = population[i].ntrlLoci[2*nNtrlLoci*j + locus*2 + 1];
			
			cont_table_tot[allele1]++;
			cont_table_tot[allele2]++;
			cont_table_demes[i][allele1]++;
			cont_table_demes[i][allele2]++;
			
			if(allele1==allele2){ // if homozygous
				Q1 += 1.0;
			}
		}
	}
	
	// Q1
	Q1 /= nTot;
	
	// Q2
	double num_Q2 = 0.0;
	double denom_Q2 = 0.0;
	double sqr_Pij = 0.0;
	for(i=0; i<nDemes; i++){
		sqr_Pij = 0.0;
		for(j=0; j<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES; j++){
			sqr_Pij += pow(cont_table_demes[i][j]/(2.0*population[i].nIndividus), 2);
		}
		num_Q2 += sqr_Pij * 2.0*population[i].nIndividus * (2.0*population[i].nIndividus-1)/2.0;
		denom_Q2 += 2.0*population[i].nIndividus * (2.0*population[i].nIndividus-1)/2.0;
	}
	Q2 = num_Q2/denom_Q2;
	
	// Q3
	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		Q3 += pow(cont_table_tot[i]/(2.0*nTot), 2);
	}
	
//	printf("Q1=%lf\nQ2=%lf\nQ3=%lf\n\n", Q1, Q2, Q3);
	
	// compute n, S1 and S2 (section 7.6 de la doc de la version 4.7 de genePop)
	// nc=(S1-S2/S1)/(n-1)
	int n = 0; // number of non-empty groups
	int S1 = 0; // total sample size
	double S2 = 0.0; // sum of squared group sizes
	
	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		if(cont_table_tot[i] > 0){
			++n;
			S1 = S1 + cont_table_tot[i];
			S2 = S2 + pow(cont_table_tot[i], 2); 
		}
	}

	
	if( n<=1 ){
		target[0] = 0.0;
	}else{
		target[0] = (S1-S2/S1)/(n-1.0);
	}
	
	target[1] = Q1;
	target[2] = Q2;
	target[3] = Q3;
	
	free(cont_table_tot);
}


void nc(const Deme* population, const int nDemes, const int nNtrlLoci, const int locus, double* target){
	// fills the array target[nc, Hs, Htot]
	int i = 0;
	int j = 0;
	long allele1 = 0;
	long allele2 = 0;

	double* cont_table_tot = NULL;
	cont_table_tot = malloc((MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1)*sizeof(double));
	
	if(cont_table_tot == NULL){
		exit(0);
	}

	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		cont_table_tot[i] = 0;
	}	

	for(i=0; i<nDemes; i++){
		for(j=0; j<population[i].nIndividus; j++){
			allele1 = population[i].ntrlLoci[2*nNtrlLoci*j + locus*2];
			allele2 = population[i].ntrlLoci[2*nNtrlLoci*j + locus*2 + 1];
			
			cont_table_tot[allele1]++;
			cont_table_tot[allele2]++;
		}
	}

	// compute Ho and Hs
	double nInd = 0.0; // number of individuals in the deme_i
	double Ho = 0.0;
	double Ho_deme_i = 0.0;
	double Hs = 0.0;
	for(i=0; i<nDemes; i++){
		nInd = 0.0;
		Ho_deme_i = 0.0;
		double* cont_table_pop = NULL;
		cont_table_pop = malloc((MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1)*sizeof(double));
		
		if(cont_table_pop == NULL){
			exit(0);
		}
		
		for(j=0; j<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; j++){
			cont_table_pop[j]=0;
		}


		for(j=0; j<population[i].nIndividus; j++){
			nInd++;
			allele1 = population[i].ntrlLoci[2*nNtrlLoci*j + locus*2];
			allele2 = population[i].ntrlLoci[2*nNtrlLoci*j + locus*2 + 1];
			if(allele1 != allele2){
				Ho_deme_i++;
			}
				
			cont_table_pop[allele1]++;
			cont_table_pop[allele2]++;
		}
		Hs += heteroZ(cont_table_pop);
		free(cont_table_pop);
		Ho += Ho_deme_i/nInd;
	}
	Hs /= nDemes; // mean Hs
	Ho /= nDemes; // observed heterozygoty in the metapop

	// compute n, S1 and S2 (section 7.6 de la doc de la version 4.7 de genePop)
	int n = 0;
	int S1 = 0;
	double S2 = 0.0;
	
	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		if(cont_table_tot[i] > 0){
			++n;
			S1 = S1 + cont_table_tot[i];
			S2 = S2 + pow(cont_table_tot[i], 2); 
		}
	}

	double Htot = 0.0;
	Htot = heteroZ(cont_table_tot);
	
	free(cont_table_tot);
	
	if( n<=1 ){
		target[0] = 0.0;
	}else{
		target[0] = (S1-S2/S1)/(n-1.0);
	}
	target[1] = Hs; // expected heterozygoty averaged over demes = 1 - sum( p_i^2 )
	target[2] = Htot; // expected heterozygoy in the whole metapopulation
	target[3] = Ho; // observed proportion of heterozygotes
}


double heteroZ(const double* cont_table){
	int i = 0;
	double n_ind = 0.0;
	double Hs = 1.0;
	
	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		n_ind += 1.0 * cont_table[i];
	} 
	
	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		Hs -= pow(cont_table[i]/n_ind, 2);
	} 

	return(Hs);	
	
}

void global_stat(Deme* population, const int nDemes, const long nNtrlLoci, double* diff_stats, double* polym_stats){
	double* nc_Q1_Q2_Q3 = NULL;
	double* nc_Hs_Htot_Ho = NULL;
	nc_Q1_Q2_Q3 = malloc(4 * sizeof(double)); // contains [nc, Q1, Q2, Q3]
	nc_Hs_Htot_Ho = malloc(4 * sizeof(double)); // contains [nc, Hs, Htot, Ho]
	
	if( nc_Q1_Q2_Q3 == NULL || nc_Hs_Htot_Ho == NULL ){
		exit(0);
	}

	int k = 0;
	int j = 0;
	int i = 0;
	int cnt = 0;
	double z_bar_j = 0.0;	
	
	double* z_bar = NULL;
	z_bar = malloc(nNtrlLoci * sizeof(double));
	if(z_bar == NULL){
		exit(0);
	}

	double* var_tot = NULL;
	var_tot = malloc(nNtrlLoci * sizeof(double));
	if(var_tot == NULL){
		exit(0);
	}

	double* var_among_patches = NULL;
	var_among_patches = malloc(nNtrlLoci * sizeof(double));
	if(var_among_patches == NULL){
		exit(0);
	}

	// global statistics
	// locus specific values of Q1, Q2, Q3
	double Q1 = 0.0;
	double Q2 = 0.0;
	double Q3 = 0.0;
	
	// locus specific values of nc, Hs, Htot and Ho
	double nc_k = 0.0;
	double Hs = 0.0;
	double Htot = 0.0;
	double Ho = 0.0;

	// global fst Charles
	double num_fst_cm = 0.0;
	double denom_fst_cm = 0.0;

	// global fst coal
	double num_fst_coal = 0.0;
	double denom_fst_coal = 0.0;
	
	// global g'st
	double num_gpst = 0.0;
	double denom_gpst = 0.0;

	// global Jost's D
	double num_D = 0.0;
	double denom_D = 0.0;

	// global Fis
	double num_Fis = 0.0;
	double denom_Fis = 0.0;
	
	// global fst ANOVA
	double num_fst_anova = 0.0;
	double denom_fst_anova = 0.0;
	
	// global fis ANOVA
	double num_fis_anova = 0.0;
	double denom_fis_anova = 0.0;
	
	// global fit ANOVA
	double num_fit_anova = 0.0;
	double denom_fit_anova = 0.0;
	
	// z_bar
	for(k=0; k<nNtrlLoci; k++){ // loop over loci k
		z_bar[k] = 0.0;
		cnt = 0;
		for(j=0; j<nDemes; j++){ // loop over demes j
			for(i=0; i<population[j].nIndividus; i++){ // loop over individuals i
				z_bar[k] += population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 0]; // allele 1
				z_bar[k] += population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 1]; // allele 2
				cnt += 2;
			}
		}
		z_bar[k] = z_bar[k]/cnt;
	}
	// variance among patches
	for(k=0; k<nNtrlLoci; k++){ // loop over loci k
		var_among_patches[k] = 0;

		cnt = 0;
		
		// compute z_bar_j
		for(j=0; j<nDemes; j++){ // loop over demes j
			z_bar_j = 0.0;
			cnt = 0;
			for(i=0; i<population[j].nIndividus; i++){ // loop over individuals i
				z_bar_j += population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 0]; // allele 1
				z_bar_j += population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 1]; // allele 2
				cnt += 2;
			}
			z_bar_j /= cnt;
			
			var_among_patches[k] += pow(z_bar_j - z_bar[k], 2);
		}
		var_among_patches[k] /= nDemes;
		
	}

	// total variance in the population
	for(k=0; k<nNtrlLoci; k++){ // loop over loci k
		nc_k = 0.0;
		Hs = 0.0; 
		Htot = 0.0;
		Ho = 0.0;
		for(j=0; j<3; j++){
			nc_Q1_Q2_Q3[j] = 0.0;
			nc_Hs_Htot_Ho[j] = 0.0;
		}

		var_tot[k] = 0.0;
		cnt = 0;
		for(j=0; j<nDemes; j++){ // loop over demes j
			for(i=0; i<population[j].nIndividus; i++){ // loop over individuals i
				var_tot[k] += pow(population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 0] - z_bar[k], 2); // allele 1
				var_tot[k] += pow(population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 1] - z_bar[k], 2); // allele 2
				cnt += 2;
			}
		}
		var_tot[k] = var_tot[k]/cnt;
		
		// with Hs Htot Ho
		nc(population, nDemes, nNtrlLoci, k, nc_Hs_Htot_Ho);
		
		// with Q1 Q2 Q3
		QoneQtwoQthree(population, nDemes, nNtrlLoci, k, nc_Q1_Q2_Q3);
		Q1 = nc_Q1_Q2_Q3[1];
		Q2 = nc_Q1_Q2_Q3[2];
		Q3 = nc_Q1_Q2_Q3[3];

		// Fst ANOVA	
		num_fst_anova += (Q2-Q3) * nc_Q1_Q2_Q3[0];
		denom_fst_anova += (1-Q3) * nc_Q1_Q2_Q3[0];
		
		// Fis ANOVA	
		num_fis_anova += (Q1-Q2) * nc_Q1_Q2_Q3[0];
		denom_fis_anova += (1-Q2) * nc_Q1_Q2_Q3[0];
		
		// Fit ANOVA	
		num_fit_anova += (Q1-Q3) * nc_Q1_Q2_Q3[0];
		denom_fit_anova += (1-Q3) * nc_Q1_Q2_Q3[0];
		
		// continue the computations	
		nc_k = nc_Hs_Htot_Ho[0];
		Hs = nc_Hs_Htot_Ho[1];
		Htot = nc_Hs_Htot_Ho[2];
		Ho = nc_Hs_Htot_Ho[3];

		// global nc
		nc_Hs_Htot_Ho[0] += nc_k;
		
		// global Hs
		nc_Hs_Htot_Ho[1] += Hs;
		
		// global Htot
		nc_Hs_Htot_Ho[2] += Htot;
		
		// global Ho
		nc_Hs_Htot_Ho[3] += Ho;
		
		// global Fst charles
		num_fst_cm += var_among_patches[k] * nc_k;
		denom_fst_cm += var_tot[k] * nc_k;
		
		// global Fst coal
		num_fst_coal += (Htot - Hs) * nc_k;
		denom_fst_coal += Htot * nc_k;
	
		// global g'st
		if( Htot != 0 ){
			num_gpst += ((Htot - Hs)/Htot) * nc_k; 
			denom_gpst += ((nDemes-1)*(1-Hs)/(nDemes-1+Hs)) * nc_k;
		}

		// global Jost's D
		num_D += (Htot - Hs) * nDemes * nc_k;
		denom_D += (1 - Hs) * (nDemes-1) * nc_k;
		
		// global Fis
		if( Hs != 0 ){
			num_Fis += (Hs - Ho) * nc_k;
			denom_Fis += Hs * nc_k;
		}
	}

	double global_fst_cm = 0.0;
	if( denom_fst_cm == 0.0){
		global_fst_cm = -9;
	}else{
		global_fst_cm = num_fst_cm / denom_fst_cm;
	}
	diff_stats[0] = global_fst_cm;

	double global_fst_coal = 0.0;
	if( denom_fst_coal == 0.0){
		global_fst_coal = -9;
	}else{
		global_fst_coal = num_fst_coal / denom_fst_coal;
	}
	diff_stats[1] = global_fst_coal;

	double global_gpst = 0.0;
	if( denom_gpst == 0.0){
		global_gpst = -9;
	}else{
		global_gpst = num_gpst / denom_gpst;
	}
	diff_stats[2] = global_gpst;

	double global_D = 0.0;
	if( denom_D == 0.0){
		global_D = -9;
	}else{
		global_D = num_D / denom_D;
	}
	diff_stats[3] = global_D;

	double global_Fis = 0.0;
	if( denom_Fis == 0.0){
		global_Fis = -9;
	}else{
		global_Fis = num_Fis / denom_Fis;
	}
	diff_stats[4] = global_Fis;

	double global_Fst_anova = 0.0;
	if( denom_fst_anova == 0.0){
		global_Fst_anova = -9;
	}else{
		global_Fst_anova = num_fst_anova / denom_fst_anova;
	}
	diff_stats[5] = global_Fst_anova; 
	
	double global_Fis_anova = 0.0;
	if( denom_fis_anova == 0.0){
		global_Fis_anova = -9;
	}else{
		global_Fis_anova = num_fis_anova / denom_fis_anova;
	}
	diff_stats[6] = global_Fis_anova; 
	
	double global_Fit_anova = 0.0;
	if( denom_fit_anova == 0.0){
		global_Fit_anova = -9;
	}else{
		global_Fit_anova = num_fit_anova / denom_fit_anova;
	}
	diff_stats[7] = global_Fit_anova; 
	
	// polymorphism
	// global nc
	polym_stats[0] += nc_k / (1.0*nNtrlLoci);
	
	// global Hs
	polym_stats[1] += Hs / (1.0*nNtrlLoci);
	
	// global Htot
	polym_stats[2] += Htot / (1.0*nNtrlLoci);
	
	// global Ho
	polym_stats[3] += Ho / (1.0*nNtrlLoci);

	// free memory
	//free(fst);
	free(z_bar);
	free(var_tot);
	free(var_among_patches);
	free(nc_Q1_Q2_Q3);
	free(nc_Hs_Htot_Ho);
}

int valueInArray(const float val, const int sizeArr, const int* array){
	// test whether the value val is in the array arr
	int i;
	for(i = 0; i < sizeArr; i++){
		if(array[i] == val){
			return 1;
		}
	}
	return 0;
}


void get_parent_for_dispersal(gsl_rng *r, const int recolonization, double deme_target_mother[], double deme_target_father[], int ind_target_mother[], int ind_target_father[], const Deme* population, const int nDemes, const int extinctionStatus[], const double selfingRate, const double dispersal, const int local_deme, const int dispersal_model){
	// entries:
	//	r
	//	number of parents to sample
	//	the deme's IDs where the sampled father come from
	//	the deme's IDs where the sampled mother come from
	//	the individual's IDs within the sampled deme of the sampled father 
	//	the individual's IDs within the sampled deme of the sampled mother 
	//	the parental metapopulation
	//	number of demes in the metapopulation
	//	extinction status
	//	selfing rate
	//	pollen dispersal rate
	//	local deme where the dispersed seed will arrived
	//	model of dispersal; 0 = migrant pool model 1 = propagule pool model
	
	// get the total number of possible fathers
	int nIndividusTotal = 0; // number of individuals in the metapop (minus extincted demes and the local deme)
	int nNonExtinctedDemes = 0; // number of non extincted demes
	int deme_tmp = 0;
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==0){ // if the deme is not extincted
			if(deme_tmp!=local_deme){ // excludes individuals from the local deme
				nIndividusTotal += population[deme_tmp].nIndividus;
				nNonExtinctedDemes += 1;
			}
		}
	}
	
	// get informations about all possible parents
	double* deme_IDs = NULL;
	double* deme_fertilities = NULL;
	
	int* deme_ind = NULL;
	int* IDind_demes = NULL;
	double* IDind_metapop = NULL;
	double* maleAlloc_metapop = NULL;
	double* femaleAlloc_metapop = NULL;


	deme_IDs = malloc(nNonExtinctedDemes * sizeof(double));	
	deme_fertilities = malloc(nNonExtinctedDemes * sizeof(double));	
	
	deme_ind = malloc(nIndividusTotal * sizeof(int)); // [deme 0; deme 0; deme 0; deme 1; deme 1; deme 1]
	IDind_demes = malloc(nIndividusTotal * sizeof(int)); // [0; 1; 2; 0; 1; 2]
	IDind_metapop = malloc(nIndividusTotal * sizeof(double)); // [0; 1; 2; 3; 4; 5]
	maleAlloc_metapop = malloc(nIndividusTotal * sizeof(double)); // [0.6; 0.7; 0.5; 0.5; 0.6; 0.9]
	femaleAlloc_metapop = malloc(nIndividusTotal * sizeof(double)); // [0.4; 0.3; 0.5; 0.5; 0.4; 0.1]

	if(deme_IDs == NULL || deme_fertilities == NULL || deme_ind == NULL || IDind_demes == NULL || IDind_metapop == NULL || maleAlloc_metapop == NULL || femaleAlloc_metapop == NULL){
		exit(0);
	}
	
	int cnt_deme = 0;
	int cnt_ind = 0;
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==0){
			if(deme_tmp!=local_deme){ // excludes individuals from the local deme
				deme_IDs[cnt_deme] = deme_tmp;
				deme_fertilities[cnt_deme] = population[deme_tmp].fertility;
				cnt_deme += 1;
				
				int ind_tmp = 0;
				for(ind_tmp=0; ind_tmp<population[deme_tmp].nIndividus; ind_tmp++){
					deme_ind[cnt_ind] = deme_tmp; // deme's id for each male
					IDind_demes[cnt_ind] = ind_tmp; // male's id within its deme
					IDind_metapop[cnt_ind] = (double)cnt_ind; // male's id within the metapopulation
					maleAlloc_metapop[cnt_ind] = population[deme_tmp].maleAllocation[ind_tmp]; // male's allocation
					femaleAlloc_metapop[cnt_ind] = population[deme_tmp].femaleAllocation[ind_tmp]; // female's allocation
					cnt_ind += 1;
				}
			}
		}
	}
	
	// MOTHERS
	// sample the demes for the k mothers
	weightedSample(r, deme_IDs, deme_fertilities, deme_target_mother, nNonExtinctedDemes, recolonization);
	if(recolonization>1){
		if(dispersal_model==1){ // 0 = migrant pool model 1 = propagule pool model
			int ind_tmp = 1;
			for(ind_tmp=1; ind_tmp<recolonization; ind_tmp++){
				deme_target_mother[ind_tmp] = deme_target_mother[0]; // First returned value 
			}
		}
	}
	// sample the k mothers within the k demes
	int seed_tmp = 0;
	for(seed_tmp=0; seed_tmp<recolonization; seed_tmp++){
		double* list_females = NULL;
		double* female_fitness = NULL;
		double* sampledMother = NULL;
		
		list_females = malloc( population[ (int)deme_target_mother[seed_tmp]].nIndividus * sizeof(double));
		female_fitness = malloc( population[ (int)deme_target_mother[seed_tmp]].nIndividus * sizeof(double));
		sampledMother = malloc( 1 * sizeof(double));
		
		int ind_tmp=0;
		for(ind_tmp=0; ind_tmp<population[ (int)deme_target_mother[seed_tmp]].nIndividus; ind_tmp++){
			list_females[ind_tmp] = ind_tmp; 
			female_fitness[ind_tmp] = population[ (int)deme_target_mother[seed_tmp]].femaleAllocation[ind_tmp];
		}
		
		weightedSample(r, list_females, female_fitness, sampledMother, population[ (int)deme_target_mother[seed_tmp]].nIndividus, 1);
		ind_target_mother[seed_tmp] = (int)sampledMother[0]; // Second returned value
		
		free(list_females);
		free(female_fitness);
		free(sampledMother);
	}
	
	// FATHERS
	// test for selfing rate and dispersal when producing colonizers
	for(seed_tmp=0; seed_tmp<recolonization; seed_tmp++){ // loop over the k migrant/colonizing individuals to produce
		int autofec_test = 0;
		
		autofec_test = gsl_ran_binomial(r, selfingRate, 1);
		if(autofec_test==1){
			// selfing
			deme_target_father[seed_tmp] = deme_target_mother[seed_tmp]; // Third returned value (if selfing)
			ind_target_father[seed_tmp] = (int)ind_target_mother[seed_tmp]; // Fourth returned value (if selfing)
		}else{
			// pollen competition
			// deme_ind = malloc(nIndividusTotal * sizeof(int)); // [deme 0; deme 0; deme 0; deme 1; deme 1; deme 1]
			// IDind_demes = malloc(nIndividusTotal * sizeof(int)); // [0; 1; 2; 0; 1; 2]
			// IDind_metapop = malloc(nIndividusTotal * sizeof(double)); // [0; 1; 2; 3; 4; 5]
			// maleAlloc_metapop = malloc(nIndividusTotal * sizeof(double)); // [0.6; 0.7; 0.5; 0.5; 0.6; 0.9]
			double* male_fitness = NULL;
			male_fitness = malloc(nIndividusTotal * sizeof(double));
			int ind_tmp = 0;
			for(ind_tmp=0; ind_tmp<nIndividusTotal; ind_tmp++){
				deme_tmp = 0;
				deme_tmp = deme_ind[ind_tmp];
				if(deme_tmp == local_deme){
					male_fitness[ind_tmp] = 0;
				}else{
					if(deme_tmp==deme_target_mother[seed_tmp]){ // if the male comes from the same deme than the previously sampled mother
						male_fitness[ind_tmp] = maleAlloc_metapop[ind_tmp] * (1-dispersal);
					}else{ // if the male comes from a deme which is different from the previously sampled mother
						male_fitness[ind_tmp] = maleAlloc_metapop[ind_tmp] * dispersal / (nDemes-1);
					}
				}
			} // end of loop over all individuals of the metapopulation
			
			double* sampledFather = NULL;
			sampledFather = malloc( 1 * sizeof(double));
			
			weightedSample(r, IDind_metapop, male_fitness, sampledFather, nIndividusTotal, 1);
			
			deme_target_father[seed_tmp] = deme_ind[ (int)sampledFather[0] ]; // Third returned value (if not selfing)
			ind_target_father[seed_tmp] = (int)IDind_demes[ (int)sampledFather[0] ]; // Fourth returned value (if not selfing)
	
			free(sampledFather);
			free(male_fitness);
		} // end of condition on self fertilization
	} // end of loop over the k fathers

	free(deme_IDs);
	free(deme_fertilities);
	
	free(deme_ind);
	free(IDind_demes);
	free(IDind_metapop);
	free(maleAlloc_metapop);
	free(femaleAlloc_metapop);
}

void replacement(Deme* population, const Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci){
/*	int nIndividus;
	long* ntrlLoci;	// contains the neutral alleles for nIndividus (number of individuals within the deme) x 2 (diploids) x nNtlrLoci (number of neutral loci)
	double* quantiLoci; // contains the allelic effects for nIndividus (number of individuals within the deme) x 2 (diploids) x nQuantiLoci (number of quantitative loci)
	double* femaleAllocation; // =sum of the allelic effects over the quantitative loci
	double* maleAllocation; // =(1 - femaleAllocation)
	double fertility; // sum of individual female allocations within the deme
*/
	int deme_tmp = 0;
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
/*		int population[deme_tmp].nIndividus = 0;
		long* population[deme_tmp].ntrlLoci = NULL;
		double* population[deme_tmp].quantiLoci = NULL;
		double* population[deme_tmp].femaleAllocation = NULL;
		double* population[deme_tmp].maleAllocation = NULL;
*/		
		population[deme_tmp].nIndividus = newPopulation[deme_tmp].nIndividus;
		population[deme_tmp].ntrlLoci = malloc(2 * newPopulation[deme_tmp].nIndividus * nNtrlLoci * sizeof(long));
		population[deme_tmp].quantiLoci = malloc(2 * newPopulation[deme_tmp].nIndividus * nQuantiLoci * sizeof(double));
		population[deme_tmp].femaleAllocation = malloc(newPopulation[deme_tmp].nIndividus * sizeof(double));
		population[deme_tmp].maleAllocation = malloc(newPopulation[deme_tmp].nIndividus * sizeof(double));

		if(population[deme_tmp].ntrlLoci == NULL || population[deme_tmp].quantiLoci == NULL || population[deme_tmp].femaleAllocation == NULL || population[deme_tmp].maleAllocation == NULL){
			exit(0);
		}

		// copy paste newPopulation into population of non extincted demes
		// // neutral alleles
		int ind_tmp = 0;
		int locus = 0;
		for(ind_tmp=0; ind_tmp<population[deme_tmp].nIndividus; ind_tmp++){
			for(locus = 0; locus<nNtrlLoci; locus++){
				int allele = 0;
				for( allele=0; allele<2; allele++){
					population[deme_tmp].ntrlLoci[(2*nNtrlLoci*ind_tmp) + (2*locus) + allele] = newPopulation[deme_tmp].ntrlLoci[(2*nNtrlLoci*ind_tmp) + (2*locus) + allele];
				}
			}
		}
		
		// // quantitative alleles
		population[deme_tmp].fertility = 0.0;
		for(ind_tmp=0; ind_tmp<population[deme_tmp].nIndividus; ind_tmp++){
			population[deme_tmp].femaleAllocation[ind_tmp] = 0.0;
			for(locus = 0; locus<nQuantiLoci; locus++){
				int allele = 0;
				for( allele=0; allele<2; allele++){
					population[deme_tmp].quantiLoci[(2*nQuantiLoci*ind_tmp) + (2*locus) + allele] = newPopulation[deme_tmp].quantiLoci[(2*nQuantiLoci*ind_tmp) + (2*locus) + allele];
					population[deme_tmp].femaleAllocation[ind_tmp] += newPopulation[deme_tmp].quantiLoci[(2*nQuantiLoci*ind_tmp) + (2*locus) + allele];
					population[deme_tmp].fertility += newPopulation[deme_tmp].quantiLoci[(2*nQuantiLoci*ind_tmp) + (2*locus) + allele];
				}
			}
			population[deme_tmp].maleAllocation[ind_tmp] = 1.0-population[deme_tmp].femaleAllocation[ind_tmp];
		}
	}
}

void afficherPopulation(const Deme* population, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int extinctionStatus[], const int generation){
	// called to print some informations about population in a debug mode
	printf("GENERATION (afficherPopulation) %d\n", generation);
	int deme_tmp = 0;
	int test_debug = 0;	
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		int ind_tmp = 0;
		for(ind_tmp=0; ind_tmp<population[deme_tmp].nIndividus; ind_tmp++){
			printf("Deme: %d Ind: %d Ntrl: ", deme_tmp, ind_tmp);
			int allele = 0;
			for(allele=0; allele<(2*nNtrlLoci); allele++){
				printf("%ld ", population[deme_tmp].ntrlLoci[2*ind_tmp*nNtrlLoci+allele]);
			}
			printf(" Quanti: ");
			allele = 0;
			for(allele=0; allele<(2*nQuantiLoci); allele++){
				printf("%.4lf ", population[deme_tmp].quantiLoci[2*ind_tmp*nQuantiLoci+allele]);
				if(population[deme_tmp].quantiLoci[2*ind_tmp*nQuantiLoci+allele]==0){
					test_debug=1;
				}
			}
		printf("femAlloc: %.4lf malAlloc: %.4lf\n", population[deme_tmp].femaleAllocation[ind_tmp], population[deme_tmp].maleAllocation[ind_tmp]);
		}
		printf("Deme: %d Fertility: %.4lf\n\n", deme_tmp, population[deme_tmp].fertility);
	}

	printf("Will be extincted at the next generation: ");	
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		if(extinctionStatus[deme_tmp]==1){
			printf("D%d\t", deme_tmp);
		}
	}
	printf("\n");
	
	if(test_debug==1){
		printf("an allele 0.0000 appeared\n");
		exit(1);
	}
}

void mutation_step(gsl_rng* r, Deme* population, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const double ntrlMutation, const double quantiMutation){
	int deme_tmp = 0;
	const double epsilon = 0.00001;
	for(deme_tmp=0; deme_tmp<nDemes; deme_tmp++){
		// mutation
		// // neutral alleles
		int allele_tmp = 0;
		for(allele_tmp=0; allele_tmp<2*nNtrlLoci*population[deme_tmp].nIndividus; allele_tmp++){
			int test_mutation = 0;
			test_mutation = gsl_ran_binomial(r, ntrlMutation, 1);
			if(test_mutation==1){
				population[deme_tmp].ntrlLoci[allele_tmp] = gsl_rng_uniform_int(r, MAX_NUMBER_OF_INITIAL_NTRL_ALLELES-1) + 1; // ntrl alleles in [1, MAX_NUMBER_OF_INITIAL_NTRL_ALLELES[
			}
		}
		
		// // quantitative alleles
		int ind_tmp = 0;
		for(ind_tmp=0; ind_tmp<population[deme_tmp].nIndividus; ind_tmp++){
			int locus = 0;
			for(locus=0; locus<nQuantiLoci; locus++){
				allele_tmp = 0;
				for(allele_tmp=0; allele_tmp<2; allele_tmp++){
					int test_mutation = 0;
					test_mutation = gsl_ran_binomial(r, quantiMutation, 1);
					if(test_mutation==1){
						double currentAllelicEffect = 0.0;
						currentAllelicEffect = population[deme_tmp].quantiLoci[ (2*nQuantiLoci*ind_tmp) + (2*locus) + allele_tmp ];
						
						double minEffect = epsilon; // minimum value that a quantitative mutation can take
						minEffect = (1 - RANGE) * currentAllelicEffect;
						if(minEffect<epsilon){
							minEffect = epsilon;
						}
						
						double maxEffect = 1.0/2.0/nQuantiLoci;	// maximum value that a quantitative mutation can take
						maxEffect = (1 + RANGE) * currentAllelicEffect;
						if( maxEffect > 1.0/2.0/nQuantiLoci){
							maxEffect = 1.0/2.0/nQuantiLoci;
						}
						
						double newAllelicEffect = 0.0;
						newAllelicEffect = gsl_ran_flat(r, minEffect, maxEffect);
						if(newAllelicEffect<epsilon){
							newAllelicEffect=epsilon;
						}
						
						// print for debug
//						printf("D:%d\tI:%d\tL:%d\tA:%d\told:%lf\tnew:%lf\n", deme_tmp, ind_tmp, locus, allele_tmp, currentAllelicEffect, newAllelicEffect);

						population[deme_tmp].quantiLoci[ (2*nQuantiLoci*ind_tmp) + (2*locus) + allele_tmp ] = newAllelicEffect;
						
						population[deme_tmp].femaleAllocation[ind_tmp] = population[deme_tmp].femaleAllocation[ind_tmp] - currentAllelicEffect + newAllelicEffect;
						population[deme_tmp].maleAllocation[ind_tmp] = 1-population[deme_tmp].femaleAllocation[ind_tmp];
						
						population[deme_tmp].fertility = population[deme_tmp].fertility - currentAllelicEffect + newAllelicEffect;
					} // end of condition on mutation
				} // end of loop over alleles 0 and 1
			} // end of loop over loci
		} // end of loop over individuals
	} // end of loop over demes
} // end of function

double z(const int maxIndPerDem, const double selfingRate, const double migrationRate, const double dispersal, const double extinction, const int recolonization){
	int n=maxIndPerDem;
	double a=selfingRate;
	double m=migrationRate/(1.0*n);
	double mp=dispersal;
	double e=extinction;
	int k=recolonization;
	
	// forward
//	return((k*(-(-1+pow(a,2))*mp*(4+(-1+a)*mp)*(-1+n)+2*m*(2+(-1+a)*mp)*((-1+pow(a,2))*mp*(-1+n)+2*(1-a+n+a*n))-pow(m,2)*(2+(-1+a)*mp)*((-1+pow(a,2))*mp*(-1+n)+2*(1-a+n+a*n)))+4*pow(e,4)*mp*((1+a)*mp*(-1+n)-2*(1+a)*m*mp*(-1+n)+pow(m,2)*((-1+a)*(-1+n)+mp*(-2+(2+k+a*k)*n)))+e*(-4*(1+a-2*m*(2+(-1+a)*mp)+pow(m,2)*(2+(-1+a)*mp))*(-1+n)+k*((-1-3*a+pow(a,2)+3*pow(a,3))*pow(mp,2)*(-1+n)+4*(1-a+n+a*n)+4*(-1+a)*mp*(-1-2*a+2*(1+a)*n)+pow(m,2)*(8*(1+a)*mp*(1+a*(-1+n))+4*(1-a+n+a*n)+(-1+a)*pow(mp,2)*(3-4*a-3*pow(a,2)+n+4*a*n+3*pow(a,2)*n))-2*m*(8*(1-a+n+a*n)+4*mp*(1+a+2*pow(a,2)*(-1+n)+2*a*n)+(-1+a)*pow(mp,2)*(3-4*a-3*pow(a,2)+n+4*a*n+3*pow(a,2)*n))))-pow(e,2)*(8-8*n+(-1+a)*k*pow(mp,2)*(-1-8*a-3*pow(a,2)+(5+8*a+3*pow(a,2))*n)+4*mp*(-(1+3*a)*(-1+n)+k*(2-a-pow(a,2)+pow((1+a),2)*n))-2*m*(4*(-3+a)*(-1+n)+4*mp*(-(1+3*a)*(-1+n)+(3+a)*k*(1-a+n+a*n))+pow(mp,2)*(4*(-1+a+n-a*n)+(-1+3*a)*k*(3-2*a-pow(a,2)+pow((1+a),2)*n)))+pow(m,2)*(-4*(2+2*(1+a)*mp+(-1+a)*pow(mp,2))*(-1+n)+k*(-4*(1-a+n+a*n)+4*mp*(3-2*a-pow(a,2)+pow((1+a),2)*n)+pow(mp,2)*(1+7*a-5*pow(a,2)-3*pow(a,3)+pow((1+a),2)*(-1+3*a)*n))))+pow(e,3)*(mp*(-4*(3+a+2*a*mp)*(-1+n)+(1+a)*k*mp*(3-2*a-pow(a,2)+pow((1+a),2)*n))+pow(m,2)*(4*(-1+a+n-a*n)-4*mp*(-4+k-a*k+4*n+2*k*n+2*a*k*n)+pow(mp,2)*(-8*a*(-1+n)+k*(7-3*a-3*pow(a,2)-pow(a,3)+pow((1+a),3)*n)))-2*m*mp*(-16*(-1+n)+mp*(4*(-1+3*a+n-3*a*n)+(1+a)*k*(3-2*a-pow(a,2)+(5+2*a+pow(a,2))*n)))))/(2*(-k*(-2*m*(2+(-1+a)*mp)*((-1+a)*mp*(-1+n)+2*n)+pow(m,2)*(2+(-1+a)*mp)*((-1+a)*mp*(-1+n)+2*n)+(-1+a)*mp*(4*n+mp*(-1-a+(-1+a)*n)))+2*pow(e,4)*mp*((1+a)*mp*(-1+n)-2*(1+a)*m*mp*(-1+n)+pow(m,2)*((-1+a)*(-1+n)+2*mp*(-1+n+k*n)))+e*(4+(-1+a)*(2+k+3*a*k)*pow(mp,2)*(-1+n)-4*n+4*k*n+2*(-1+a)*mp*(2-2*n+k*(-1+4*n))-2*m*(4+(-4+8*k)*n+(-1+a)*k*pow(mp,2)*(3-3*a+n+3*a*n)+2*mp*(-1+a+n-a*n+k*(3-3*a+4*a*n)))+pow(m,2)*(4+4*(-1+k)*n+(-1+a)*k*pow(mp,2)*(1-3*a+n+3*a*n)+2*mp*(-1+a+n-a*n+k*(2-2*a+4*a*n))))-pow(e,2)*(4-4*n+(-1+a)*k*pow(mp,2)*(-1-3*a+(5+3*a)*n)+2*mp*(-(1+3*a)*(-1+n)+k*(1-a+2*(1+a)*n))+pow(m,2)*(4*(1+a)*mp*(1+(-1+k)*n)-4*(-1+n+k*n)+pow(mp,2)*(2*(-1+a+n-a*n)+(1+a)*k*(3-3*a-n+3*a*n)))-2*m*(4-4*n+pow(mp,2)*(4*(-1+a+n-a*n)+(-1+3*a)*k*(1-a+n+a*n))+2*mp*(-(3+a)*(-1+n)+2*k*(1-a+(3+a)*n))))+pow(e,3)*mp*(-2*(3+a+2*a*mp)*(-1+n)+(1+a)*k*mp*(1-a+n+a*n)+pow(m,2)*(-2*(1+a)*(2+mp)*(-1+n)+k*(-2+2*a-8*n+mp*(3-2*a-pow(a,2)+pow((1+a),2)*n)))-2*m*(8-8*n+mp*(-2*(-1+3*a)*(-1+n)+k*(1-pow(a,2)+(5+2*a+pow(a,2))*n)))))));

	// backward
	return((8*pow(e,2)* (-1+m)*(-1+n)+(1+a)*k*pow(m,2)* pow((2+(-1+a)*mp),2)* (-1+n)-4* e* (-1-a+2* m+(-1+a)* mp)* (-1+n)+(-1+a)* k* mp* ((-1+pow(a,2))* mp* (-1+n)+4* (-2-a+n+a* n))-2* k* m* (4* (-1+pow(a,2))* mp* (-1+n)+pow((-1+a),2)* (1+a)* pow(mp,2)* (-1+n)+4*(-a+n+a*n))-e*k*(-1+m)*(-4*(-1+pow(a,2))*mp* (-1+n)-pow((-1+a),2)* (1+a)* pow(mp,2)* (-1+n)+(1+a)* m* pow((2+(-1+a)* mp),2) *(-1+n)-4 *(1-a+n+a* n)))/(8 *pow(e,2)* (-1+m)* (-1+n)+2* k* (pow(m,2)* pow((2+(-1+a)* mp),2)* (-1+n)+(-1+a)* mp* (4+(-1+a)* mp)* (-1+n)-2* m* (-2+4* (-1+a)* mp* (-1+n)+pow((-1+a),2)*pow(mp,2)* (-1+n)+4* n))-2* e* (-1+m)* (4* (-1+n)+k* (-pow((-1+a),2)* pow(mp,2)* (-1+n)+m* pow((2+(-1+a)* mp),2)* (-1+n)-4* n+4* mp* (-1+a+n-a* n)))));
}

