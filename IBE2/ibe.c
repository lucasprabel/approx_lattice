#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "arithmetic.h"
#include "sampling.h"



/*
	Generates the master public key A and the master secret key (R, cplx_R, sch_comp).
*/
void Setup(poly_matrix A, poly_matrix R, cplx_poly_matrix cplx_R, cplx_poly_matrix sch_comp)
	{

	scalar *r_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly_matrix r = r_coeffs;

	scalar *e_coeffs = malloc(PARAM_N * PARAM_M * sizeof(scalar));
	poly_matrix e = e_coeffs;

		
	// Sampling of r and e
	SampleR_centered((signed_poly_matrix) r, PARAM_SIGMA);

	SampleR_matrix_centered((signed_poly_matrix) e, PARAM_M, 1, PARAM_SIGMA);

	// Compute the trapdoor R

	zero_poly(R, PARAM_N * (PARAM_M+1) * PARAM_M - 1);

	scalar *e_add_inverse_coeffs = malloc(PARAM_N * PARAM_M * sizeof(scalar));
	poly_matrix e_add_inverse = e_add_inverse_coeffs;
	for (int i = 0 ; i < PARAM_N * PARAM_M ; i++){
		e_add_inverse[i] = -e[i];
	}


	memcpy(R, e_add_inverse, PARAM_N * PARAM_M * sizeof(scalar));
	for (int i = 1 ; i <= PARAM_M ; i++){
		memcpy(&R[i * PARAM_M * PARAM_N  + (i-1) * PARAM_N], r, PARAM_N * sizeof(scalar));
	}

	// Compute the Schur complements (that are necessary to the Gaussian preimage sampling operation) + cplx_R
	construct_complex_private_key(cplx_R, sch_comp, R);
	
	// Add q to each component of R (so that R's coeffs are positive) and put it in the CRT domain
	for(int i = 0 ; i < PARAM_N * (PARAM_M+1) * PARAM_M ; i++)
		{
		R[i] += PARAM_Q;
		}

	matrix_crt_representation(R, PARAM_M + 1, PARAM_M, LOG_R);

	//Compute r^{-1}

	scalar *r_inv_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly_matrix r_inv = r_inv_coeffs;
	invert_poly(r_inv, r, PARAM_N, 1);
	crt_representation(r_inv, LOG_R);

	//Compute f+e

	scalar *f_coeffs = malloc((PARAM_K-PARAM_L) * PARAM_N * sizeof(scalar));
	poly_matrix f = f_coeffs;

	for(int i = 0; i < PARAM_K - PARAM_L ; i++){
		poly f_i = poly_matrix_element(f, 1, i, 0);
		f_i[0] = pow(PARAM_B, PARAM_L + i);
	}


	// Compute a = [1 | a']
	add_to_poly_matrix(e, f, PARAM_M, 1);

	scalar* a_prime_coeffs = malloc(2*PARAM_N * PARAM_M * sizeof(scalar));
	poly_matrix a_prime = a_prime_coeffs;
	zero_poly(a_prime, 2*PARAM_N * PARAM_M-1);

	for (int i = 0 ; i < PARAM_M ; i ++){
		double_poly a_prime_i = (double_poly) poly_matrix_element(a_prime, 1, i, 0);
		poly e_i = poly_matrix_element(e, 1, i, 0);
		mul_crt_poly(a_prime_i, r_inv, e_i, LOG_R);
	}

	zero_poly(A, PARAM_N * (PARAM_M + 1) - 1);

	poly A_0 = poly_matrix_element(A, PARAM_M + 1, 0, 0);
	A_0[0] = 1;

	memcpy(A + PARAM_N, a_prime, PARAM_N * PARAM_M * sizeof(scalar));

	//Reduce A's coefficients mod q
	freeze_poly(A, PARAM_N * (PARAM_M + 1)-1);

	
	free(r_coeffs);
	free(e_coeffs);
	//free(a_prime);
	}


/*
	Generate the secret key x of an identity id using the master secret key (R, cplx_R, sch_comp) and the master public key (A,u).
*/
void Extract(poly_matrix x, poly_matrix A, poly_matrix R, cplx_poly_matrix cplx_R, cplx_poly_matrix sch_comp, scalar *id)
	{
	// Compute id's inverse and put it in the CRT domain
	scalar id_inv_coeffs[PARAM_N];
	poly id_inv = id_inv_coeffs;
	
	invert_poly(id_inv, id, PARAM_N, 1);
	crt_representation(id_inv, LOG_R);
	
	// Use id to construct A_id
	//construct_A_m(A, id);
	
	// Sample x
	// Use of approx_sample_pre_target, which is approximate sample_pre with a target u (needs id_inv as an argument)
	approx_sample_pre_2(x, A, R, cplx_R, sch_comp, id_inv);
	
	// Deconstruct A_m
	//deconstruct_A_m(A, id);
	coeffs_representation(id, LOG_R);
	}
	

void Encrypt(poly_matrix A, scalar *id, poly M, poly_matrix b, poly c)
{

	// Compute id's inverse and put it in the CRT domain
	scalar id_inv_coeffs[PARAM_N];
	poly id_inv = id_inv_coeffs;
	
	invert_poly(id_inv, id, PARAM_N, 1);
	crt_representation(id_inv, LOG_R);

	// Variables
	scalar *s_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly_matrix s = s_coeffs;

	scalar *e_1_coeffs = malloc(PARAM_M * PARAM_N * sizeof(scalar));
	scalar *e_2_coeffs = malloc(PARAM_N * sizeof(scalar));

	poly_matrix e_1 = e_1_coeffs, e_2 = e_2_coeffs;

	// Sampling of s
	SampleR_centered((signed_poly_matrix) s, PARAM_TAU);

	// Sampling of the errors

	//e_1 <- D_{R^{m,1},tau}
	SampleR_matrix_centered((signed_poly_matrix) e_1, PARAM_M, 1, PARAM_TAU);

	// e_2 <- D_{R,tau}
	SampleR_centered((signed_poly_matrix) e_2, PARAM_TAU);

	

	// Compute b
	memset(b, 0, PARAM_M * PARAM_N * sizeof(scalar));

	//crt_representation(s, LOG_R);
	//mul_crt_poly_matrix(b, s, A, 1, 1, PARAM_M, LOG_R);
	matrix_coeffs_representation(b, PARAM_M, 1, LOG_R);

	//add_to_poly_matrix(b, e_1, PARAM_M, 1);

	// Compute c
	memset(c, 0, PARAM_N * sizeof(scalar));

	//mul_crt_poly((double_poly) c, id_inv, s, LOG_R);
	coeffs_representation(c, LOG_R);
	//add_poly(c, c, e_2, PARAM_N-1);	

	for (int i=0 ; i < PARAM_N ; i++)
		{
			M[i] = floor(PARAM_Q /2) * M[i];
		}

	add_poly(c, c, M, PARAM_N-1);

	free(s);
	free(e_1);
	free(e_2);
}



void Decrypt(poly_matrix x, poly_matrix b, poly c, poly M)
{
	scalar *res_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly res = res_coeffs;

	scalar *factor_coeffs = malloc(PARAM_N * sizeof(scalar));
	poly_matrix factor = factor_coeffs;

	matrix_crt_representation(b, PARAM_M, 1, LOG_R);

	for (int i = 0 ; i < PARAM_N * (PARAM_M+1) ; i++)
	{
		x[i] += PARAM_Q;
	}

	matrix_crt_representation(x, PARAM_M+1, 1, LOG_R);
	mul_crt_poly_matrix(factor, b, &x[PARAM_N], 1, PARAM_M, 1, LOG_R);
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
			M[i] = 1;
		}

		else {
			M[i] = 0;
		}
	}

}