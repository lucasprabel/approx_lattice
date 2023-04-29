#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "arithmetic.h"
#include "sampling.h"
#include "cpucycles.h"

extern unsigned long long timing_sampleZ_Setup;
extern unsigned long long timing_sampleZ_Extract;
extern unsigned long long timing_sample_s;
extern unsigned long long timing_sample_e0;
extern unsigned long long timing_sample_e1;
extern unsigned long long timing_sample_e_prime;



/*
	Generates the master public key (A,u) and the master secret key (T, cplx_T, sch_comp).
*/
void Setup(poly_matrix A, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, poly_matrix u)
	{

	unsigned long long begin_timing1 = 0;
	unsigned long long end_timing1 = 0;
	// scalar A_hat_coeffs[PARAM_D * PARAM_D * PARAM_N], AprimeT_coeffs[PARAM_D * PARAM_D * PARAM_K * PARAM_N];

	scalar *A_barre_coeffs = malloc(PARAM_D * PARAM_M_BARRE * PARAM_N * sizeof(scalar)), *AbarreT_coeffs = malloc(PARAM_D * PARAM_W * PARAM_N * sizeof(scalar));
	poly_matrix A_barre = A_barre_coeffs, AbarreT = AbarreT_coeffs;
	
	// A_barre <- U(R_q^{d,bar(m)}) is considered to be in the CRT domain
	random_poly(A_barre, PARAM_D * PARAM_M_BARRE * PARAM_N - 1);
		
	// T <- D_{R^{2d,d(k-l)},sigma}
	begin_timing1 = cpucycles_start();
	SampleR_matrix_centered((signed_poly_matrix) T, PARAM_M_BARRE, PARAM_W, PARAM_SIGMA);
	end_timing1 = cpucycles_stop();
	timing_sampleZ_Setup += (end_timing1 - begin_timing1);

	// Compute the Schur complements (that are necessary to the Gaussian preimage sampling operation) + cplx_T
	construct_complex_private_key(cplx_T, sch_comp, T);
	
	// Add q to each component of T (so that T's coeffs are positive) and put it in the CRT domain
	for(int i = 0 ; i < PARAM_N * PARAM_M_BARRE * PARAM_W ; ++i)
		{
		T[i] += PARAM_Q;
		}

	matrix_crt_representation(T, PARAM_M_BARRE, PARAM_W, LOG_R);

	mul_crt_poly_matrix(AbarreT, A_barre, T, PARAM_D, PARAM_M_BARRE, PARAM_W, LOG_R);
	
	// A = (A_barre | -A_barre*T)
	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_i0 = poly_matrix_element(A, PARAM_M, i, 0);
		poly_matrix A_barre_i = poly_matrix_element(A_barre, PARAM_M_BARRE, i, 0);
		
		memcpy(A_i0, A_barre_i, PARAM_M_BARRE * PARAM_N * sizeof(scalar));
		}

	for(int i = 0 ; i < PARAM_D ; ++i)
		{
		poly_matrix A_i1 = poly_matrix_element(A, PARAM_M, i, PARAM_M_BARRE);
		poly_matrix AbarreT_i = poly_matrix_element(AbarreT, PARAM_W, i, 0);

		for(int j = 0 ; j < PARAM_W * PARAM_N ; ++j)
			{
			A_i1[j] = 2*PARAM_Q - AbarreT_i[j];
			}
		}
	
	// Reduce A's coefficients mod q
	freeze_poly(A, PARAM_N * PARAM_D * PARAM_M - 1);


	random_poly(u, PARAM_N * PARAM_D - 1);

	
	free(A_barre);
	free(AbarreT);
	}



/*
	Generate the secret key of an identity id using the master secret key (T, cplx_T, sch_comp) and the master public key (A,u).
*/
void Extract(poly_matrix x, poly_matrix A, poly_matrix u, poly_matrix T, cplx_poly_matrix cplx_T, cplx_poly_matrix sch_comp, scalar *id)
	{
	unsigned long long begin_timing1 = 0;
	unsigned long long end_timing1 = 0;
	
	// Compute id's inverse and put it in the CRT domain
	scalar id_inv_coeffs[PARAM_N];
	poly id_inv = id_inv_coeffs;
	
	invert_poly(id_inv, id, PARAM_N, 1);
	crt_representation(id_inv, LOG_R);
	
	// Use id to construct A_id
	construct_A_m(A, id);
	// Sample x
	// Use of approx_sample_pre_target, which is approximate sample_pre with a target u (needs id_inv as an argument)
	begin_timing1 = cpucycles_start();
	approx_sample_pre_target(x, A, T, cplx_T, sch_comp, id_inv, u);
	end_timing1 = cpucycles_stop();
	timing_sampleZ_Extract += (end_timing1 - begin_timing1);

	// Deconstruct A_m
	deconstruct_A_m(A, id);
	}
/*
	
*/
void Encrypt(poly_matrix A, poly_matrix u, scalar *id, poly M, poly_matrix b, poly c)
{
	unsigned long long begin_timing1 = 0;
	unsigned long long end_timing1 = 0;
	unsigned long long begin_timing2 = 0;
	unsigned long long end_timing2 = 0;
	unsigned long long begin_timing3 = 0;
	unsigned long long end_timing3 = 0;
	unsigned long long begin_timing4 = 0;
	unsigned long long end_timing4 = 0;




	// Use id to construct A_id
	construct_A_m(A, id);


	scalar *s_coeffs = malloc(PARAM_D * PARAM_N * sizeof(scalar));
	poly_matrix s = s_coeffs;

	scalar *e_0_coeffs = malloc(PARAM_M_BARRE * PARAM_N * sizeof(scalar));
	scalar *e_1_coeffs = malloc(PARAM_W * PARAM_N * sizeof(scalar));
	scalar *e_prime_coeffs = malloc(PARAM_N * sizeof(scalar));

	poly_matrix e_0 = e_0_coeffs, e_1 = e_1_coeffs;
	poly e_prime = e_prime_coeffs;


	begin_timing1 = cpucycles_start();
	// Sampling of s
	SampleR_matrix_centered((signed_poly_matrix) s, PARAM_D, 1, PARAM_TAU);
	end_timing1 = cpucycles_stop();
	timing_sample_s += (end_timing1 - begin_timing1);

	// Sampling of the errors

	// e_0 <- D_{R^{m'-k,1},tau}
	begin_timing2 = cpucycles_start();
	SampleR_matrix_centered((signed_poly_matrix) e_0, PARAM_M_BARRE, 1, PARAM_TAU);
	end_timing2 = cpucycles_stop();
	timing_sample_e0 += (end_timing2 - begin_timing2);

	// e_1 <- D_{R^{k,1},gamma}
	begin_timing3 = cpucycles_start();
	SampleR_matrix_centered((signed_poly_matrix) e_1, PARAM_W, 1, PARAM_GAMMA);
	end_timing3 = cpucycles_stop();
	timing_sample_e1 += (end_timing3 - begin_timing3);

	// e_prime <- D_{R,tau}
	begin_timing4 = cpucycles_start();
	SampleR_centered((signed_poly_matrix) e_prime, PARAM_TAU);
	end_timing4 = cpucycles_stop();
	timing_sample_e_prime += (end_timing4 - begin_timing4);




	memset(b, 0, PARAM_M * PARAM_N * sizeof(scalar));

	// Compute s^T * A_id
	scalar *factor1_coeffs = malloc(PARAM_M * PARAM_N * sizeof(scalar)); // (s^T * A_id)
	poly_matrix factor1 = factor1_coeffs;

	// Make sure s has positive coefficients and put it in the CRT domain
	for(int i = 0 ; i < PARAM_N * PARAM_D; ++i)
	{
		s[i] += PARAM_Q;
	}
	
	matrix_crt_representation(s, PARAM_D, 1, LOG_R);
	multiply_by_A_transpose_approx(factor1, A, s);
	matrix_coeffs_representation(factor1, PARAM_M, 1, LOG_R);


	// Compute (e_0^T | e_1^T)^T
	scalar *factor2_coeffs = malloc(PARAM_M * PARAM_N * sizeof(scalar));
	poly_matrix factor2 = factor2_coeffs;
	memcpy(factor2, e_0, PARAM_M_BARRE * PARAM_N * sizeof(scalar));
	memcpy(factor2 + PARAM_M_BARRE * PARAM_N , e_1, PARAM_W * PARAM_N * sizeof(scalar));

	add_to_poly_matrix(b, factor1, 1, PARAM_M);
	add_to_poly_matrix(b, factor2, 1, PARAM_M);

	//matrix_crt_representation(u, PARAM_D, 1, LOG_R);
	mul_crt_poly_matrix(c, s, u, 1, PARAM_D, 1, LOG_R);
	coeffs_representation(c, LOG_R);

	add_poly(c, c, e_prime, PARAM_N-1);

	for (int i=0 ; i < PARAM_N ; i++)
		{
			M[i] = floor(PARAM_Q /2) * M[i];
		}

	add_poly(c, c, M, PARAM_N-1);
	freeze_poly(c, PARAM_N-1);

	deconstruct_A_m(A, id);

	free(s);
	free(e_0);
	free(e_1);
	free(e_prime);
}



void Decrypt(poly_matrix x, poly_matrix b, poly c, poly M)
{
	scalar *res_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly res = res_coeffs;

	scalar *factor_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly_matrix factor = factor_coeffs;

	// Make sure b has positive coefficients and put it in the CRT domain
	for(int i = 0 ; i < PARAM_M * PARAM_N; ++i)
	{
		b[i] += PARAM_Q;
	}
	matrix_crt_representation(b, PARAM_M, 1, LOG_R);

	for (int i = 0 ; i < PARAM_N * PARAM_M ; i++) // useful ?
	{
		x[i] += PARAM_Q;
	}

	// matrix_crt_representation(x, PARAM_M, 1, LOG_R); x deja sous forme CRT
	mul_crt_poly_matrix(factor, b, x, 1, PARAM_M, 1, LOG_R);
	matrix_coeffs_representation(factor, 1, 1, LOG_R);

	for (int i = 0 ; i < PARAM_N ; i++)
	{
		factor[i] = 2*PARAM_Q + (c[i] - factor[i]);
	}

	freeze_poly(factor, PARAM_N-1);

	res = factor;


	for (int i = 0 ; i < PARAM_N ; i++)
	{
		if (res[i] < floor(PARAM_Q/2)){
			M[i] = 0;
		}

		else {
			M[i] = 1;
		}
	}

}

