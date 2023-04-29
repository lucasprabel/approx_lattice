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

#define MESSAGE_BYTES 512 // keep it divisible by 16
#define NTESTS 100
#define CPU_CYCLES (4.2 * 1000000000.0)

unsigned long long timing_overhead;
unsigned long long timing_sampleZ_G = 0;
unsigned long long timing_sampleZ_P = 0;
unsigned long long timing_sampleG = 0;
unsigned long long timing_samplePerturb = 0;
unsigned long long timing_sampleArith = 0;

unsigned long long timing_extract = 0;
unsigned long long timing_encrypt = 0;
unsigned long long timing_decrypt = 0;

unsigned long long timing_sampleZ_Setup = 0;
unsigned long long timing_precomp_Setup = 0;
unsigned long long timing_arith_Setup = 0;
unsigned long long timing_setup = 0;

void time_setup(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * (PARAM_M + 1) * sizeof(scalar)), *R_coeffs = malloc(PARAM_N * (PARAM_M + 1) * PARAM_M * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * (PARAM_M + 1)* (PARAM_M+2) / 2 * sizeof(cplx)), *cplx_R_coeffs = malloc(PARAM_N * (PARAM_M + 1) * PARAM_M * sizeof(cplx));
	poly_matrix A = A_coeffs, R = R_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_R = cplx_R_coeffs;



	for(unsigned i = 0; i < NTESTS; ++i) {

		begin_timing = cpucycles_start();
		
		// Generate Keys
		Setup(A, R, cplx_R, sch_comp);
		
		end_timing = cpucycles_stop();
		timing_setup += (end_timing - begin_timing);
 	}

 	printf("----------- Setup -----------\n");
 	timing_setup = timing_setup/NTESTS - timing_overhead;
 	printf("Total: %lld cycles (%.2lf ms)\n", timing_setup, (timing_setup*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(R);
	free(sch_comp);
	free(cplx_R);
	}

void time_extract(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	scalar *A_coeffs = malloc(PARAM_N * (PARAM_M + 1) * sizeof(scalar)), *R_coeffs = malloc(PARAM_N * (PARAM_M + 1) * PARAM_M * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * (PARAM_M + 1)* (PARAM_M+2) / 2 * sizeof(cplx)), *cplx_R_coeffs = malloc(PARAM_N * (PARAM_M + 1) * PARAM_M * sizeof(cplx));
	poly_matrix A = A_coeffs, R = R_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_R = cplx_R_coeffs;

	// Generate Keys
	Setup(A, R, cplx_R, sch_comp);



	for(unsigned i = 0; i < NTESTS; ++i) {
		// Generate an identity
		scalar id_coeffs[PARAM_N] = {0};
		poly id = id_coeffs;
		
		random_poly(id, SMALL_DEGREE - 1);

 		// Compute a signature (nu)
		scalar x_coeffs[PARAM_N * (PARAM_M + 1)];
		poly_matrix x = x_coeffs;

		begin_timing = cpucycles_start();
		Extract(x, A, R, cplx_R, sch_comp, id);
		end_timing = cpucycles_stop();
		timing_extract += (end_timing - begin_timing);
 	}

 	printf("\n----------- Extract -----------\n");
 	timing_extract = timing_extract/NTESTS - timing_overhead;
 	printf("Extract: %lld cycles (%.2lf ms)\n", timing_extract, (timing_extract*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(R);
	free(sch_comp);
	free(cplx_R);
	}

void time_encrypt(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * (PARAM_M + 1) * sizeof(scalar)), *R_coeffs = malloc(PARAM_N * (PARAM_M + 1) * PARAM_M * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * (PARAM_M + 1)* (PARAM_M+2) / 2 * sizeof(cplx)), *cplx_R_coeffs = malloc(PARAM_N * (PARAM_M + 1) * PARAM_M * sizeof(cplx));
	poly_matrix A = A_coeffs, R = R_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_R = cplx_R_coeffs;

	// Generate Keys
	Setup(A, R, cplx_R, sch_comp);
	
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
		Encrypt(A, id, m, b, c);
		end_timing = cpucycles_stop();
		timing_encrypt += (end_timing - begin_timing);
 	}

 	printf("\n----------- Encrypt -----------\n");
 	timing_encrypt = timing_encrypt/NTESTS - timing_overhead;
 	printf("Encrypt: %lld cycles (%.2lf ms)\n", timing_encrypt, (timing_encrypt*1000)/CPU_CYCLES);
 	
 	printf("\n\n");

 	free(A);
	free(R);
	free(sch_comp);
	free(cplx_R);
	}

void time_decrypt(void)
	{
	timing_overhead = cpucycles_overhead();
	unsigned long long begin_timing = 0;
	unsigned long long end_timing = 0;

	
	scalar *A_coeffs = malloc(PARAM_N * (PARAM_M + 1) * sizeof(scalar)), *R_coeffs = malloc(PARAM_N * (PARAM_M + 1) * PARAM_M * sizeof(scalar));
	cplx *sch_comp_coeffs = malloc(PARAM_N * (PARAM_M + 1)* (PARAM_M+2) / 2 * sizeof(cplx)), *cplx_R_coeffs = malloc(PARAM_N * (PARAM_M + 1) * PARAM_M * sizeof(cplx));
	poly_matrix A = A_coeffs, R = R_coeffs;
	cplx_poly_matrix sch_comp = sch_comp_coeffs, cplx_R = cplx_R_coeffs;

	// Generate Keys
	Setup(A, R, cplx_R, sch_comp);

	// Generate a message
	scalar m_coeffs[PARAM_N] = {0};
	poly m = m_coeffs;
	
	random_message(m, PARAM_N-1);

	// Generate an identity
	scalar id_coeffs[PARAM_N] = {0};
	poly id = id_coeffs;
	
	random_poly(id, SMALL_DEGREE - 1);

	// Generate nu = sk_id
	scalar nu_coeffs[PARAM_N * (PARAM_M+1)];
	poly_matrix nu = nu_coeffs;
	random_poly(nu, PARAM_N * (PARAM_M+1)-1);


	//Extract(nu, A, u, T, cplx_T, sch_comp, id);

 	// Compute b and c
 	scalar *b_coeffs = malloc(PARAM_N * PARAM_M * sizeof(scalar));
 	scalar *c_coeffs = malloc(PARAM_N * sizeof(scalar));

 	poly_matrix b = b_coeffs;
 	poly c = c_coeffs;


	
	Encrypt(A, id, m, b, c);


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
	free(R);
	free(sch_comp);
	free(cplx_R);
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
	printf("q           : %d\n", PARAM_Q);
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

