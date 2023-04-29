#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "random.h"
#include "sampling.h"
#include "ibe.h"
#include "arithmetic.h"

#include "cpucycles.h"

#define NTESTS 1 // If NTESTS is too low, timings will be higher
#define CPU_CYCLES (4.2 * 1000000000.0)

unsigned long long timing_overhead;
unsigned long long timing_sampleZ_G = 0;
unsigned long long timing_sampleZ_P = 0;
unsigned long long timing_sampleG = 0;
unsigned long long timing_samplePerturb = 0;
unsigned long long timing_sampleArith = 0;

unsigned long long timing_decrypt = 0;

unsigned long long timing_setup = 0;
unsigned long long timing_sampleZ_Setup = 0;

unsigned long long timing_extract = 0;
unsigned long long timing_sampleZ_Extract = 0;

unsigned long long timing_encrypt = 0;
unsigned long long timing_sampleZ_Encrypt = 0;
unsigned long long timing_sample_s = 0;
unsigned long long timing_sample_e0 = 0;
unsigned long long timing_sample_e1 = 0;
unsigned long long timing_sample_e_prime = 0;

void time_setup(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * PARAM_M * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * PARAM_M_BARRE * PARAM_W * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_M_BARRE * (PARAM_M_BARRE + 1) / 2 * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * PARAM_M_BARRE * PARAM_W * sizeof(cplx));
	scalar *u_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	poly_matrix A = A_coeffs, T = T_coeffs, u = u_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;



	for(unsigned i = 0; i < NTESTS; ++i) {

		begin_timing = cpucycles_start();
		
		// Generate Keys
		Setup(A, T, cplx_T, sch_comp, u);
		
		end_timing = cpucycles_stop();
		timing_setup += (end_timing - begin_timing);
 	}

 	timing_sampleZ_Setup = timing_sampleZ_Setup/NTESTS - timing_overhead;
 	timing_setup = timing_setup/NTESTS - timing_overhead;


 	printf("----------- Setup -----------\n");
 	printf("Gaussian Sampling      : %lld cycles (%.2lf ms) (%.2lf %% of Setup)\n", timing_sampleZ_Setup, (timing_sampleZ_Setup*1000)/CPU_CYCLES, ((double) timing_sampleZ_Setup/ (double) timing_setup)*100.0);
 	printf("Total: %lld cycles (%.2lf ms)\n", timing_setup, (timing_setup*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	free(u);
	}

void time_extract(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * PARAM_M * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * PARAM_M_BARRE * PARAM_W * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_M_BARRE * (PARAM_M_BARRE + 1) / 2 * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * PARAM_M_BARRE * PARAM_W * sizeof(cplx));
	scalar *u_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	poly_matrix A = A_coeffs, T = T_coeffs, u = u_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;

	// Generate Keys
	Setup(A, T, cplx_T, sch_comp, u);



	for(unsigned i = 0; i < NTESTS; ++i) {
		// Generate an identity
		scalar id_coeffs[PARAM_N] = {0};
		poly id = id_coeffs;
		
		random_poly(id, SMALL_DEGREE - 1);

 		// Compute a signature (nu)
		scalar nu_coeffs[PARAM_N * PARAM_M];
		poly_matrix nu = nu_coeffs;

		begin_timing = cpucycles_start();
		Extract(nu, A, u, T, cplx_T, sch_comp, id);
		end_timing = cpucycles_stop();
		timing_extract += (end_timing - begin_timing);
 	}

	timing_sampleZ_Extract = timing_sampleZ_Extract/NTESTS - timing_overhead;
	timing_extract = timing_extract/NTESTS - timing_overhead;


 	printf("\n----------- Extract -----------\n");
	printf("Preimage Sampling     : %lld cycles (%.2lf ms) (%.2lf %% of Extract)\n", timing_sampleZ_Extract, (timing_sampleZ_Extract*1000)/CPU_CYCLES, ((double) timing_sampleZ_Extract/ (double) timing_extract)*100.0);
 	printf("Extract: %lld cycles (%.2lf ms)\n", timing_extract, (timing_extract*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	free(u);
	}

void time_encrypt(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * PARAM_M * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * PARAM_M_BARRE * PARAM_W * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_M_BARRE * (PARAM_M_BARRE + 1) / 2 * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * PARAM_M_BARRE * PARAM_W * sizeof(cplx));
	scalar *u_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	poly_matrix A = A_coeffs, T = T_coeffs, u = u_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;

	// Generate Keys
	Setup(A, T, cplx_T, sch_comp, u);
	
	for(unsigned i = 0; i < NTESTS; ++i) {

		// Generate a message
		scalar m_coeffs[PARAM_N] = {0};
		poly m = m_coeffs;
		
		random_message(m, PARAM_N-1);

		// Generate an identity
		scalar id_coeffs[PARAM_N] = {0};
		poly id = id_coeffs;
	
		random_poly(id, SMALL_DEGREE - 1);

 		// Compute b and c
 		scalar *b_coeffs = malloc(PARAM_N * PARAM_M * sizeof(scalar));
 		scalar *c_coeffs = malloc(PARAM_N * sizeof(scalar));

 		poly_matrix b = b_coeffs;
 		poly c = c_coeffs;


		begin_timing = cpucycles_start();
		Encrypt(A, u, id, m, b, c);
		end_timing = cpucycles_stop();
		timing_encrypt += (end_timing - begin_timing);
 	}

	timing_sample_s = timing_sample_s/NTESTS - timing_overhead;
	timing_sample_e0 = timing_sample_e0/NTESTS - timing_overhead;
	timing_sample_e1 = timing_sample_e1/NTESTS - timing_overhead;
	timing_sample_e_prime = timing_sample_e_prime/NTESTS - timing_overhead;
	timing_encrypt = timing_encrypt/NTESTS - timing_overhead;
	timing_sampleZ_Encrypt = timing_sample_s + timing_sample_e0 + timing_sample_e1 + timing_sample_e_prime;
 	printf("\n----------- Encrypt -----------\n");
	printf("Gaussian Sampling of s        : %lld cycles (%.2lf ms) (%.2lf %% of Encrypt)\n", timing_sample_s, (timing_sample_s*1000)/CPU_CYCLES, ((double) timing_sample_s/ (double) timing_encrypt)*100.0);
	printf("Gaussian Sampling of e0       : %lld cycles (%.2lf ms) (%.2lf %% of Encrypt)\n", timing_sample_e0, (timing_sample_e0*1000)/CPU_CYCLES, ((double) timing_sample_e0/ (double) timing_encrypt)*100.0);
	printf("Gaussian Sampling of e1       : %lld cycles (%.2lf ms) (%.2lf %% of Encrypt)\n", timing_sample_e1, (timing_sample_e1*1000)/CPU_CYCLES, ((double) timing_sample_e1/ (double) timing_encrypt)*100.0);
	printf("Gaussian Sampling of e_prime  : %lld cycles (%.2lf ms) (%.2lf %% of Encrypt)\n", timing_sample_e_prime, (timing_sample_e_prime*1000)/CPU_CYCLES, ((double) timing_sample_e_prime/ (double) timing_encrypt)*100.0);
	printf("Total Gaussian Sampling       : %lld cycles (%.2lf ms) (%.2lf %% of Encrypt)\n", timing_sampleZ_Encrypt, (timing_sampleZ_Encrypt*1000)/CPU_CYCLES, ((double) timing_sampleZ_Encrypt/ (double) timing_encrypt)*100.0);
 	printf("Encrypt: %lld cycles (%.2lf ms)\n", timing_encrypt, (timing_encrypt*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	free(u);
	}

void time_decrypt(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * PARAM_D * PARAM_M * sizeof(scalar)), *T_coeffs = malloc(PARAM_N * PARAM_M_BARRE * PARAM_W * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * PARAM_M_BARRE * (PARAM_M_BARRE + 1) / 2 * sizeof(cplx)), *cplx_T_coeffs = malloc(PARAM_N * PARAM_M_BARRE * PARAM_W * sizeof(cplx));
	scalar *u_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	poly_matrix A = A_coeffs, T = T_coeffs, u = u_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_T = cplx_T_coeffs;

	// Generate Keys
	Setup(A, T, cplx_T, sch_comp, u);

	// Generate a message
	scalar m_coeffs[PARAM_N] = {0};
	poly m = m_coeffs;
	
	random_message(m, PARAM_N-1);

	// Generate an identity
	scalar id_coeffs[PARAM_N] = {0};
	poly id = id_coeffs;
	
	random_poly(id, SMALL_DEGREE - 1);

	// Generate nu = sk_id
	scalar nu_coeffs[PARAM_N * PARAM_M];
	poly_matrix nu = nu_coeffs;

	Extract(nu, A, u, T, cplx_T, sch_comp, id);

 	// Compute b and c
 	scalar *b_coeffs = malloc(PARAM_N * PARAM_M * sizeof(scalar));
 	scalar *c_coeffs = malloc(PARAM_N * sizeof(scalar));

 	poly_matrix b = b_coeffs;
 	poly c = c_coeffs;
	
	Encrypt(A, u, id, m, b, c);


	scalar *M_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly M = M_coeffs;

	for(unsigned i = 0; i < NTESTS; ++i) {


		begin_timing = cpucycles_start();
		
		Decrypt(nu, b, c, M);
		
		end_timing = cpucycles_stop();
		timing_decrypt += (end_timing - begin_timing);
 	}


 	printf("----------- Decrypt -----------\n");
 	
 	timing_decrypt = timing_decrypt/NTESTS - timing_overhead;
 	printf("Total: %lld cycles (%.2lf ms)\n\n\n", timing_decrypt, (timing_decrypt*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(T);
	free(sch_comp);
	free(cplx_T);
	free(u);
	free(b);
	free(c);
	free(M);
	}




int main(void) {
	
	init_crt_trees();
	init_D_lattice_coeffs();
	init_cplx_roots_of_unity();
	random_bytes_init();

 	printf("----------- Parameters -----------\n");
	printf("n           : %d\n", PARAM_N);
	printf("q           : %ld\n", PARAM_Q);
	printf("d           : %d\n", PARAM_D);
	printf("l           : %d\n", PARAM_L);
	printf("Sigma       : %f\n", PARAM_SIGMA);
	printf("Alpha       : %f\n", PARAM_ALPHA);
	printf("Zeta        : %f\n", PARAM_ZETA);
	printf("Tau         : %f\n", PARAM_TAU);
	printf("Gamma       : %f\n", PARAM_GAMMA);
 	printf("----------------------------------\n\n");


	
 	time_setup();
 	
 	time_extract();
 	
 	time_encrypt();

 	time_decrypt(); 

	return 0;
}

