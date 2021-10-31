/*
 *  Copyright 2021 The ZKrypt Project. All Rights Reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License"); you may
 *  not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bn254.h"
#include "bn254_test.h"


int test_bn(void)
{
	uint32_t a[8];
	uint32_t b[8];
	uint32_t r[16];

	printf("%s\n", __FUNCTION__);

	bn_from_hex(a, 8, hex_a);
	bn_from_hex(b, 8, hex_b);
	bn_print(" a = ", a, 8);
	bn_print(" b = ", b, 8);

	bn_add(r, a, b, 8);
	printf(" a + b = %s\n", hex_add_a_b);
	bn_print("       = ", r, 8);
	if (!bn_equ_hex(r, 8, hex_add_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_sub(r, a, b, 8);
	printf(" a - b = %s\n", hex_sub_a_b);
	bn_print("       = ", r, 8);
	if (!bn_equ_hex(r, 8, hex_sub_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_sub(r, b, a, 8);
	printf(" b - a = %s\n", hex_sub_b_a);
	bn_print("       = ", r, 8);
	if (!bn_equ_hex(r, 8, hex_sub_b_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_mul_lo(r, a, b, 8);
	printf(" a * b (mod 2^256) = %s\n", hex_mul_lo_a_b);
	bn_print("                   = ", r, 8);
	if (!bn_equ_hex(r, 8, hex_mul_lo_a_b)) {
		fprintf(stderr, "%s %d: error\n",  __FILE__, __LINE__);
	}

	bn_mul(r, a, b, 8);
	printf(" a * b = %s\n", hex_mul_a_b);
	bn_print("       = ", r, 16);
	if (!bn_equ_hex(r, 16, hex_mul_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	return 1;
}

int test_bn_mod(void)
{


	uint32_t a[8];
	uint32_t b[8];
	uint32_t n[8];
	uint32_t r[16];

	printf("%s\n", __FUNCTION__);

	bn_from_hex(a, 8, hex_a);
	bn_from_hex(b, 8, hex_b);
	bn_from_hex(n, 8, hex_n);

	bn_print(" a = ", a, 8);
	bn_print(" b = ", b, 8);
	bn_print(" n = ", n, 8);


	bn_mod_add_non_const_time(r, a, b, n, 8);
	bn_print(" a + b (mod n) = ", r, 8);
	bn_mod_add_const_time(r, a, b, n, 8);
	bn_print("               = ", r, 8);
	  printf("               = %s\n", hex_fn_add_a_b);
	if (!bn_equ_hex(r, 8, hex_fn_add_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_mod_sub_non_const_time(r, a, b, n, 8);
	bn_print(" a - b (mod n) = ", r, 8);
	bn_mod_sub_const_time(r, a, b, n, 8);
	bn_print("               = ", r, 8);
	  printf("               = %s\n", hex_fn_sub_a_b);
	if (!bn_equ_hex(r, 8, hex_fn_sub_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_mod_sub_non_const_time(r, b, a, n, 8);
	bn_print(" b - a (mod n) = ", r, 8);
	bn_mod_sub_const_time(r, b, a, n, 8);
	bn_print("               = ", r, 8);
	  printf("               = %s\n", hex_fn_sub_b_a);
	if (!bn_equ_hex(r, 8, hex_fn_sub_b_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_mod_neg(r, a, n, 8);
	bn_print(" -a (mod n) = ", r, 8);
	  printf("            = %s\n", hex_fn_neg_a);
	if (!bn_equ_hex(r, 8, hex_fn_neg_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	return 1;
}


int test_bn_barrett(void)
{


	uint32_t r[8];
	uint32_t a[8];
	uint32_t b[8];
	uint32_t p[8];
	uint32_t p_mu[9];
	uint32_t n[8];
	uint32_t n_mu[9];

	printf("%s\n", __FUNCTION__);

	bn_from_hex(a, 8, hex_a);
	bn_from_hex(b, 8, hex_b);
	bn_from_hex(p, 8, hex_p);
	bn_from_hex(p_mu, 9, hex_p_mu);
	bn_from_hex(n, 8, hex_n);
	bn_from_hex(n_mu, 9, hex_n_mu);

	bn_print(" a = ", a, 8);
	bn_print(" b = ", b, 8);
	bn_print(" p = ", n, 8);
	bn_print(" mu(p) = ", n_mu, 9);

	bn_barrett_mod_mul(r, a, b, p, p_mu, 8);
	bn_print(" a * b (mod p) = ", r, 8);
	printf(  "               = %s\n", hex_fp_mul_a_b);
	if (!bn_equ_hex(r, 8, hex_fp_mul_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_barrett_mod_sqr(r, a, p, p_mu, 8);
	bn_print(" a^2 (mod p) = ", r, 8);
	printf(  "             = %s\n", hex_fp_sqr_a);
	if (!bn_equ_hex(r, 8, hex_fp_sqr_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_barrett_mod_exp(r, a, b, p, p_mu, 8);
	bn_print(" a^b (mod p) = ", r, 8);
	printf(  "             = %s\n", hex_fp_exp_a_b);
	if (!bn_equ_hex(r, 8, hex_fp_exp_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_barrett_mod_inv(r, a, p, p_mu, 8);
	bn_print(" a^-1 (mod p) = ", r, 8);
	printf(  "              = %s\n", hex_fp_inv_a);
	if (!bn_equ_hex(r, 8, hex_fp_inv_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}


	bn_print(" n = ", n, 8);
	bn_print(" mu(n) = ", n_mu, 8);

	bn_barrett_mod_mul(r, a, b, n, n_mu, 8);
	bn_print(" a * b (mod n) = ", r, 8);
	printf(  "               = %s\n", hex_fn_mul_a_b);
	if (!bn_equ_hex(r, 8, hex_fn_mul_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_barrett_mod_sqr(r, a, n, n_mu, 8);
	bn_print(" a^2 (mod n) = ", r, 8);
	printf(  "             = %s\n", hex_fn_sqr_a);
	if (!bn_equ_hex(r, 8, hex_fn_sqr_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_barrett_mod_exp(r, a, b, n, n_mu, 8);
	bn_print(" a^b (mod n) = ", r, 8);
	printf(  "             = %s\n", hex_fn_exp_a_b);
	if (!bn_equ_hex(r, 8, hex_fn_exp_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_barrett_mod_inv(r, a, n, n_mu, 8);
	bn_print(" a^-1 (mod n) = ", r, 8);
	printf(  "              = %s\n", hex_fn_inv_a);
	if (!bn_equ_hex(r, 8, hex_fn_inv_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	return 1;
}

int test_bn_montgomery(void)
{
	uint32_t a[8];
	uint32_t b[8];
	uint32_t p[8];
	uint32_t p_inv_neg[8];
	uint32_t p_one_sqr[8];
	uint32_t n[8];
	uint32_t n_inv_neg[8];
	uint32_t n_one_sqr[8];
	uint32_t r[8];

	printf("%s\n", __FUNCTION__);

	bn_from_hex(a, 8, hex_a);
	bn_from_hex(b, 8, hex_b);
	bn_from_hex(p, 8, hex_p);
	bn_from_hex(p_inv_neg, 8, hex_p_inv_neg);
	bn_from_hex(p_one_sqr, 8, hex_p_one_sqr);

	bn_print(" a = ", a, 8);
	bn_print(" b = ", b, 8);
	bn_print(" p = ", p, 8);

	bn_mod_mul_montgomery(r, a, b, p_one_sqr, p, p_inv_neg, 8);
	bn_print(" a * b (mod p) = ", r, 8);
	printf(  "               = %s\n", hex_fp_mul_a_b);
	if (!bn_equ_hex(r, 8, hex_fp_mul_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_mod_sqr_montgomery(r, a, p_one_sqr, p, p_inv_neg, 8);
	bn_print(" a^2 (mod p) = ", r, 8);
	printf(  "             = %s\n", hex_fp_sqr_a);
	if (!bn_equ_hex(r, 8, hex_fp_sqr_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_from_hex(n, 8, hex_n);
	bn_from_hex(n_inv_neg, 8, hex_n_inv_neg);
	bn_from_hex(n_one_sqr, 8, hex_n_one_sqr);

	bn_print(" n = ", n, 8);

	bn_mod_mul_montgomery(r, a, b, n_one_sqr, n, n_inv_neg, 8);
	bn_print(" a * b (mod n) = ", r, 8);
	printf(  "               = %s\n", hex_fn_mul_a_b);
	if (!bn_equ_hex(r, 8, hex_fn_mul_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	bn_mod_sqr_montgomery(r, a, n_one_sqr, n, n_inv_neg, 8);
	bn_print(" a^2 (mod n) = ", r, 8);
	printf(  "             = %s\n", hex_fn_sqr_a);
	if (!bn_equ_hex(r, 8, hex_fn_sqr_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	return 1;
}



int test_fn(void)
{
	uint32_t a[8];
	uint32_t b[8];
	uint32_t r[8];
	printf("%s\n", __FUNCTION__);

	bn_from_hex(a, 8, hex_a);
	bn_from_hex(b, 8, hex_b);

	bn_print(" a = ", a, 8);
	bn_print(" b = ", b, 8);
	bn_print(" n = ", BN254_N, 8);


	fn_from_hex(a, hex_a);
	fn_from_hex(b, hex_b);

	fn_add(r, a, b);
	fn_get_bn(r, r);
	bn_print(" a + b (mod n) = ", r, 8);
	printf(  "               = %s\n", hex_fn_add_a_b);
	if (!bn_equ_hex(r, 8, hex_fn_add_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fn_sub(r, a, b);
	fn_get_bn(r, r);
	bn_print(" a - b (mod n) = ", r, 8);
	printf(  "               = %s\n", hex_fn_sub_a_b);
	if (!bn_equ_hex(r, 8, hex_fn_sub_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fn_sub(r, b, a);
	fn_get_bn(r, r);
	bn_print(" b - a (mod p) = ", r, 8);
	printf(  "               = %s\n", hex_fn_sub_b_a);
	if (!bn_equ_hex(r, 8, hex_fn_sub_b_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fn_neg(r, a);
	fn_get_bn(r, r);
	bn_print(" -a (mod p) = ", r, 8);
	printf(  "            = %s\n", hex_fn_neg_a);
	if (!bn_equ_hex(r, 8, hex_fn_neg_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fn_mul(r, a, b);
	fn_get_bn(r, r);
	bn_print(" a * b (mod n) = ", r, 8);
	printf(  "               = %s\n", hex_fn_mul_a_b);
	if (!bn_equ_hex(r, 8, hex_fn_mul_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fn_sqr(r, a);
	fn_get_bn(r, r);
	bn_print(" a^2 (mod n) = ", r, 8);
	printf(  "             = %s\n", hex_fn_sqr_a);
	if (!bn_equ_hex(r, 8, hex_fn_sqr_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fn_inv(r, a);
	fn_get_bn(r, r);
	bn_print(" a^-1 (mod n) = ", r, 8);
	printf(  "              = %s\n", hex_fn_inv_a);
	if (!bn_equ_hex(r, 8, hex_fn_inv_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	return 0;
}

int test_fp(void)
{

	uint32_t a[8];
	uint32_t b[8];
	uint32_t r[8];
	printf("%s\n", __FUNCTION__);

	bn_from_hex(a, 8, hex_a);
	bn_from_hex(b, 8, hex_b);

	bn_print(" a = ", a, 8);
	bn_print(" b = ", b, 8);

	fp_from_hex(a, hex_a);
	fp_from_hex(b, hex_b);

	fp_add(r, a, b);
	fp_get_bn(r, r);
	bn_print(" a + b (mod p) = ", r, 8);
	printf(  "               = %s\n", hex_fp_add_a_b);
	if (!bn_equ_hex(r, 8, hex_fp_add_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fp_sub(r, a, b);
	fp_get_bn(r, r);
	bn_print(" a - b (mod p) = ", r, 8);
	printf(  "               = %s\n", hex_fp_sub_a_b);
	if (!bn_equ_hex(r, 8, hex_fp_sub_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fp_sub(r, b, a);
	fp_get_bn(r, r);
	bn_print(" b - a (mod p) = ", r, 8);
	printf(  "               = %s\n", hex_fp_sub_b_a);
	if (!bn_equ_hex(r, 8, hex_fp_sub_b_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fp_neg(r, a);
	fp_get_bn(r, r);
	bn_print(" -a (mod p) = ", r, 8);
	printf(  "            = %s\n", hex_fp_neg_a);
	if (!bn_equ_hex(r, 8, hex_fp_neg_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fp_mul(r, a, b);
	fp_get_bn(r, r);
	bn_print(" a * b (mod p) = ", r, 8);
	printf(  "               = %s\n", hex_fp_mul_a_b);
	if (!bn_equ_hex(r, 8, hex_fp_mul_a_b)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fp_sqr(r, a);
	fp_get_bn(r, r);
	bn_print(" a^2 (mod p) = ", r, 8);
	printf(  "             = %s\n", hex_fp_sqr_a);
	if (!bn_equ_hex(r, 8, hex_fp_sqr_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	fp_inv(r, a);
	fp_get_bn(r, r);
	bn_print(" a^-1 (mod p) = ", r, 8);
	printf(  "              = %s\n", hex_fp_inv_a);
	if (!bn_equ_hex(r, 8, hex_fp_inv_a)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	return 0;
}

int test_point(void)
{
	point_t G1;
	point_t P;
	fn_t a;

	printf("%s\n", __FUNCTION__);

	bn_from_hex(a, 8, hex_a);
	bn_print(" a = ", a, 8);

	point_from_hex(&G1, hex_G1);
	point_print(" G1 = ", &G1);
	point_print("    = ", &BN254_G1);
	printf(     "    = %s\n", hex_G1);

	point_dbl(&P, &G1);
	point_print(" 2 * G1 = ", &P);
	printf(     "        = %s\n", hex_2G1);
	if (!point_equ_hex(&P, hex_2G1)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	point_add(&P, &P, &G1);
	point_print(" 3 * G1 = ", &P);
	printf(     "        = %s\n", hex_3G1);
	if (!point_equ_hex(&P, hex_3G1)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	point_mul_affine_non_const_time(&P, a, (affine_point_t *)&G1);
	point_print(" a * G1 = ", &P);
	if (!point_equ_hex(&P, hex_aG1)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	point_mul_affine_const_time(&P, a, (affine_point_t *)&G1);
	point_print("        = ", &P);
	if (!point_equ_hex(&P, hex_aG1)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	printf(     "        = %s\n", hex_aG1);
	point_mul_generator(&P, a);
	point_print("        = ", &P);
	if (!point_equ_hex(&P, hex_aG1)) {
		fprintf(stderr, "%s %d: error\n", __FILE__, __LINE__);
	}

	return 0;
}

int test_point_multi_mul(void)
{
	int x[8] = {1, 2, 4, 8, 16, 32, 64, 128};
	point_t P[8];
	point_t T[255];
	point_t R;

	fn_t a[8];

	printf("%s\n", __FUNCTION__);



	int i;
	for (i = 0; i < 8; i++) {
		uint32_t a[8];
		bn_set_word(a, x[i], 8);
		point_mul_generator(&P[i], a);
		//point_print(" >> ", &P[i]);
	}


	point_multi_mul_affine_pre_compute(P, T);
	/*
	for (i = 0; i < 255; i++) {
		printf("T[%d] ", i+1);
		point_print("", &T[i]);
	}
	*/
	/*
	fn_t *hex_coffs = poly_new(8);
	hex_coffs[0] = fn_from_hex(hex_fa0);
	hex_coffs[1] = fn_from_hex(hex_fa1);
	hex_coffs[2] = fn_from_hex(hex_fa2);
	hex_coffs[3] = fn_from_hex(hex_fa3);
	hex_coffs[4] = fn_from_hex(hex_fa4);
	hex_coffs[5] = fn_from_hex(hex_fa5);
	hex_coffs[6] = fn_from_hex(hex_fa6);
	hex_coffs[7] = fn_from_hex(hex_fa7);
	*/
	/*
	for (i = 0; i < 8; i++) {
		bn_from_hex(a[i], 8, hex_coffs[i]);
	}
	*/

	point_multi_mul_affine_with_pre_compute(&R, a, T);

	point_print("R = ", &R);

	return 0;
}


int test_poly(void)
{

	fn_t *fa = poly_new(16 + 4);
	fn_t *fb = poly_new(16 + 4);
	fn_t *fr = poly_new(32 + 4);
	int fa_len = 16;
	int fb_len = 16;
	int fr_len = 0;

	fn_t a;
	fn_t r;

	/*
	FILE *fp = NULL;
	if (!(fp = fopen("poly.txt", "r"))) {
		fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
		return -1;
	}
	poly_from_file(fa, 16, fp);
	poly_from_file(fb, 16, fp);
	*/

	fn_from_hex(fa[0], hex_fa0);
	fn_from_hex(fa[1], hex_fa1);
	fn_from_hex(fa[2], hex_fa2);
	fn_from_hex(fa[3], hex_fa3);
	fn_from_hex(fa[4], hex_fa4);
	fn_from_hex(fa[5], hex_fa5);
	fn_from_hex(fa[6], hex_fa6);
	fn_from_hex(fa[7], hex_fa7);
	fn_from_hex(fa[8], hex_fa8);
	fn_from_hex(fa[9], hex_fa9);
	fn_from_hex(fa[10], hex_fa10);
	fn_from_hex(fa[11], hex_fa11);
	fn_from_hex(fa[12], hex_fa12);
	fn_from_hex(fa[13], hex_fa13);
	fn_from_hex(fa[14], hex_fa14);
	fn_from_hex(fa[15], hex_fa15);

	fn_from_hex(fb[0], hex_fb0);
	fn_from_hex(fb[1], hex_fb1);
	fn_from_hex(fb[2], hex_fb2);
	fn_from_hex(fb[3], hex_fb3);
	fn_from_hex(fb[4], hex_fb4);
	fn_from_hex(fb[5], hex_fb5);
	fn_from_hex(fb[6], hex_fb6);
	fn_from_hex(fb[7], hex_fb7);
	fn_from_hex(fb[8], hex_fb8);
	fn_from_hex(fb[9], hex_fb9);
	fn_from_hex(fb[10], hex_fb10);
	fn_from_hex(fb[11], hex_fb11);
	fn_from_hex(fb[12], hex_fb12);
	fn_from_hex(fb[13], hex_fb13);
	fn_from_hex(fb[14], hex_fb14);
	fn_from_hex(fb[15], hex_fb15);


	fn_from_hex(a, hex_a);

	poly_print("fa = ", fa, 16);
	poly_print("fb = ", fb, 16);

	poly_add(fr, &fr_len, fa, fa_len, fb, fb_len);
	poly_print("fa + fb = ", fr, fr_len);

	poly_sub(fr, &fr_len, fa, fa_len, fb, fb_len);
	poly_print("fa - fb = ", fr, fr_len);

	poly_mul(fr, &fr_len, fa, fa_len, fb, fb_len);
	poly_print("fa * fb = ", fr, fr_len);

	poly_add_scalar(fr, &fr_len, fa, fa_len, a);
	poly_print("fa + a = ", fr, fr_len);

	poly_sub_scalar(fr, &fr_len, fa, fa_len, a);
	poly_print("fa - a = ", fr, fr_len);

	poly_mul_scalar(fr, &fr_len, fa, fa_len, a);
	poly_print("fa * a = ", fr, fr_len);

	poly_mul_vector(fr, &fr_len, fa, fa_len, fb, fb_len);
	poly_print("mul_vec(fa, fb) = ", fr, fr_len);


	poly_copy(fr, &fr_len, fa, fa_len);
	poly_add_blind(fr, &fr_len, fb, 2);
	poly_print("blind2(a) = ", fr, fr_len);

	poly_copy(fr, &fr_len, fa, fa_len);
	poly_add_blind(fr, &fr_len, fb, 3);
	poly_print("blind3(a) = ", fr, fr_len);

	// 不是所有的fa都可以被(x-a)
	poly_copy(fr, &fr_len, fa, fa_len);
	poly_div_x_sub_scalar(fr, &fr_len, a);
	poly_print("fa/(x-a) = ", fr, fr_len);

	// 这个计算结果似乎
	poly_copy(fr, &fr_len, fa, fa_len);
	poly_div_ZH(fr, &fr_len, fa, fa_len, 16);
	poly_print("fa/(x^n-1) = ", fr, fr_len);

	poly_eval(r, fa, fa_len, a);
	fn_print("fa(a) = ", r);

	/*

	fn_t *H = poly_new(7);
	fn_from_hex(H[0], hex_w0);
	fn_from_hex(H[1], hex_w1);
	fn_from_hex(H[2], hex_w2);
	fn_from_hex(H[3], hex_w3);
	fn_from_hex(H[4], hex_w4);
	fn_from_hex(H[5], hex_w5);
	fn_from_hex(H[6], hex_w6);

	poly_copy(fr, &fr_len, fa, 7);
	fft(fr, fr + 8, H, 1, 7);

	poly_print("fft(fa) = ", fr, 7);
	*/

	return 1;
}

int test_fft(void)
{
	fn_t *a = poly_new(16);
	fn_t *H = poly_new(16);

	fn_from_hex(a[0], hex_fa0);
	fn_from_hex(a[1], hex_fa1);
	fn_from_hex(a[2], hex_fa2);
	fn_from_hex(a[3], hex_fa3);
	fn_from_hex(a[4], hex_fa4);
	fn_from_hex(a[5], hex_fa5);
	fn_from_hex(a[6], hex_fa6);
	fn_from_hex(a[7], hex_fa7);
	fn_from_hex(a[8], hex_fa8);
	fn_from_hex(a[9], hex_fa9);
	fn_from_hex(a[10], hex_fa10);
	fn_from_hex(a[11], hex_fa11);
	fn_from_hex(a[12], hex_fa12);
	fn_from_hex(a[13], hex_fa13);
	fn_from_hex(a[14], hex_fa14);
	fn_from_hex(a[15], hex_fa15);

	fn_from_hex(H[0], hex_w0);
	fn_from_hex(H[1], hex_w1);
	fn_from_hex(H[2], hex_w2);
	fn_from_hex(H[3], hex_w3);
	fn_from_hex(H[4], hex_w4);
	fn_from_hex(H[5], hex_w5);
	fn_from_hex(H[6], hex_w6);
	fn_from_hex(H[7], hex_w7);
	fn_from_hex(H[8], hex_w8);
	fn_from_hex(H[9], hex_w9);
	fn_from_hex(H[10], hex_w10);
	fn_from_hex(H[11], hex_w11);
	fn_from_hex(H[12], hex_w12);
	fn_from_hex(H[13], hex_w13);
	fn_from_hex(H[14], hex_w14);
	fn_from_hex(H[15], hex_w15);

	fft(a, H, 4);

	poly_print("a(H) = ", a, 16);


	poly_interpolate(a, H, 4);
	poly_print("a = ", a, 16);
	return 0;
}

int main(void)
{

	test_bn();
	test_bn_mod();
	test_bn_barrett();
	test_bn_montgomery();
	test_fp();
	test_fn();
	test_point();
	test_point_multi_mul();
	test_poly();
	test_fft();

	return 0;
}
