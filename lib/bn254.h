/*
 *  Copyright 2021 The ZKrypt Project. All Rights Reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License"); you may
 *  not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 */

#ifndef ZKRYPT_BN254_H
#define ZKRYPT_BN254_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


extern const uint32_t BN254_P[8];
extern const uint32_t BN254_N[8];
extern const uint32_t BN254_P_MU[9];
extern const uint32_t BN254_N_MU[9];
extern const uint32_t BN254_P_INV_NEG[8];
extern const uint32_t BN254_P_MONT_ONE[8];
extern const uint32_t BN254_P_MONT_ONE_SQR[8];
extern const uint32_t BN254_N_INV_NEG[8];
extern const uint32_t BN254_N_MONT_ONE[8];
extern const uint32_t BN254_N_MONT_ONE_SQR[8];
extern const uint32_t BN254_ONE[8];
extern const uint32_t BN254_TWO[8];


//#define BN254_CONST_TIME_
#define BN254_FP_USE_MONTGOMERY
//#define BN254_FN_USE_MONTGOMERY

void bn_print(const char *str, const uint32_t *a, int k);
void bn_set_word(uint32_t *r, uint32_t a, int k);
#define bn_set_zero(r,k) bn_set_word(r,0,k)
#define bn_set_one(r,k) bn_set_word(r,1,k)
void bn_copy(uint32_t *r, const uint32_t *a, int k);
int bn_cmp(const uint32_t *a, const uint32_t *b, int k);
int bn_is_zero(const uint32_t *a, int k);
int bn_is_one(const uint32_t *a, int k);
int bn_add(uint32_t *r, const uint32_t *a, const uint32_t *b, int k);
int bn_sub(uint32_t *r, const uint32_t *a, const uint32_t *b, int k);
void bn_mul(uint32_t *r, const uint32_t *a, const uint32_t *b, int k);
void bn_mul_lo(uint32_t *r, const uint32_t *a, const uint32_t *b, int k);
void bn_to_bytes(const uint32_t *a, int k, uint8_t *out);
void bn_from_bytes(uint32_t *a, int k, const uint8_t *in);
void bn_to_hex(const uint32_t *a, int k, char *out);
int bn_from_hex(uint32_t *r, int k, const char *hex);
int bn_equ_hex(const uint32_t *a, int k, const char *hex);
int bn_rand(uint32_t *r, int k);
int bn_rand_range(uint32_t *r, const uint32_t *range, int k);

void bn_mod_add_non_const_time(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p, int k);
void bn_mod_sub_non_const_time(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p, int k);
void bn_mod_add_const_time(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p, int k);
void bn_mod_sub_const_time(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p, int k);
#ifdef BN254_CONST_TIME
# define bn_mod_add(r,a,b,p,k) bn_mod_add_const_time(r,a,b,p,k)
# define bn_mod_sub(r,a,b,p,k) bn_mod_sub_const_time(r,a,b,p,k)
#else
# define bn_mod_add(r,a,b,p,k) bn_mod_add_non_const_time(r,a,b,p,k)
# define bn_mod_sub(r,a,b,p,k) bn_mod_sub_non_const_time(r,a,b,p,k)
#endif
void bn_mod_neg(uint32_t *r, const uint32_t *a, const uint32_t *p, int k);

void bn_barrett_mod_mul(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p, const uint32_t *p_mu, int k);
void bn_barrett_mod_sqr(uint32_t *r, const uint32_t *a, const uint32_t *p, const uint32_t *u, int k);
void bn_barrett_mod_exp(uint32_t *r, const uint32_t *a, const uint32_t *e, const uint32_t *p, const uint32_t *p_mu, int k);
void bn_barrett_mod_inv(uint32_t *r, const uint32_t *a, const uint32_t *p, const uint32_t *p_mu, int k);

void bn_mont_mod_mul(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p, const uint32_t *p_inv_neg, int k);
void bn_mont_mod_sqr(uint32_t *r, const uint32_t *a, const uint32_t *p, const uint32_t *p_inv_neg, int k);
void bn_mont_mod_exp(uint32_t *r, const uint32_t *a, const uint32_t *e, const uint32_t *p, const uint32_t *p_inv_neg, int k);
void bn_mont_mod_inv(uint32_t *r, const uint32_t *a, const uint32_t *p, const uint32_t *p_inv_neg, int k);
void bn_mont_set(uint32_t *r, const uint32_t *a, const uint32_t *one_sqr, const uint32_t *p, const uint32_t *p_inv_neg, int k);
void bn_mont_get(uint32_t *r, const uint32_t *a, const uint32_t *p, const uint32_t *p_inv_neg, int k);

void bn_mod_mul_montgomery(uint32_t *r, const uint32_t *a, const uint32_t *b,
	const uint32_t *R_sqr, const uint32_t *p, const uint32_t *p_inv_neg, int k);
void bn_mod_sqr_montgomery(uint32_t *r, const uint32_t *a,
	const uint32_t *R_sqr, const uint32_t *p, const uint32_t *p_inv_neg, int k);


typedef uint32_t fn_t[8];

#define fn_print(s,a)		bn_print(s,a,8)
#define fn_rand(r)		bn_rand_range(r,BN254_N,8)
#define fn_copy(r,a)		bn_copy(r,a,8)
#define fn_set_word(r,a)	bn_set_word(r,a,8)
#define fn_set_zero(r)		bn_set_zero(r,8)
#define fn_is_zero(a)		bn_is_zero(a,8)
#define fn_add(r,a,b)		bn_mod_add(r,a,b,BN254_N,8)
#define fn_sub(r,a,b)		bn_mod_sub(r,a,b,BN254_N,8)
#define fn_neg(r,a)		bn_mod_neg(r,a,BN254_N,8)
#define fn_dbl(r,a)		fn_add(r,a,a)
#define fn_to_bytes(a,out)	bn_to_bytes(a,8,out)
#define fn_from_bytes(a,in)	bn_from_bytes(a,8,in)
#ifdef BN254_FN_USE_MONTGOMERY
# define fn_set_bn(r,a)		bn_mont_mod_mul(r,a,BN254_N_MONT_ONE_SQR,BN254_N,BN254_N_INV_NEG,8)
# define fn_set_one(r)		bn_copy(r,BN254_N_MONT_ONE,8)
# define fn_to_hex(a,s)		fn_get_bn(a,a); bn_to_hex(a,8,s)
# define fn_from_hex(r,s)	bn_from_hex(r,8,s); fn_set_bn(r,r)
# define fn_get_bn(a,r) 	bn_mont_get(r,a,BN254_N,BN254_N_INV_NEG,8)
# define fn_is_one(a)		(bn_cmp(a,BN254_N_MONT_ONE,8) == 0)
# define fn_mul(r,a,b)		bn_mont_mod_mul(r,a,b,BN254_N,BN254_N_INV_NEG,8)
# define fn_sqr(r,a)		bn_mont_mod_sqr(r,a,BN254_N,BN254_N_INV_NEG,8)
#else
# define fn_set_bn(r,a)		fn_copy(r,a)
# define fn_set_one(r)		bn_set_one(r,8)
# define fn_to_hex(a,s)		bn_to_hex(a,8,s)
# define fn_from_hex(r,s)	bn_from_hex(r,8,s)
# define fn_get_bn(a,r) 	fn_copy(r,a)
# define fn_is_one(a)		bn_is_one(a,8)
# define fn_mul(r,a,b)		bn_barrett_mod_mul(r,a,b,BN254_N,BN254_N_MU,8)
# define fn_sqr(r,a)		bn_barrett_mod_sqr(r,a,BN254_N,BN254_N_MU,8)
#endif

void fn_tri(fn_t r, const fn_t a);
void fn_exp(fn_t r, const fn_t a, const uint32_t *e, int elen);
void fn_inv(uint32_t *r, const uint32_t *a);
void fn_mul_word(fn_t r, const fn_t a, uint32_t word);


typedef uint32_t fp_t[8];

#define fp_print(s,a)		bn_print(s,a,8)
#define fp_copy(r,a)		bn_copy(r,a,8)
#define fp_set_zero(r)		bn_set_zero(r,8)
#define fp_rand			bn_rand_range(r,BN254_P,8)
#define fp_is_zero(a)		bn_is_zero(a,8)
#define fp_add(r,a,b)		bn_mod_add(r,a,b,BN254_P,8)
#define fp_sub(r,a,b)		bn_mod_sub(r,a,b,BN254_P,8)
#define fp_neg(r,a)		bn_mod_neg(r,a,BN254_P,8)
#define fp_dbl(r,a)		fp_add(r,a,a)
#define fp_to_bytes(a,out)	bn_to_bytes(a,8,out)
#define fp_from_bytes(a,in)	bn_from_bytes(a,8,in)
#ifdef BN254_FP_USE_MONTGOMERY
# define fp_set_bn(r,a)		bn_mont_mod_mul(r,a,BN254_P_MONT_ONE_SQR,BN254_P,BN254_P_INV_NEG,8)
# define fp_set_one(r)		bn_copy(r,BN254_P_MONT_ONE,8)
# define fp_from_hex(r,s)	bn_from_hex(r,8,s); fp_set_bn(r,r);
# define fp_get_bn(a,r) 	bn_mont_get(r,a,BN254_P,BN254_P_INV_NEG,8)
# define fp_is_one(a)		(bn_cmp(a,BN254_P_MONT_ONE,8) == 0)
# define fp_mul(r,a,b)		bn_mont_mod_mul(r,a,b,BN254_P,BN254_P_INV_NEG,8)
# define fp_sqr(r,a)		bn_mont_mod_sqr(r,a,BN254_P,BN254_P_INV_NEG,8)
#else
# define fp_set_bn(r,a)		fp_copy(r,a)
# define fp_set_one(r)		bn_set_one(r,8)
# define fp_from_hex(r,s)	bn_from_hex(r,8,s)
# define fp_get_bn(a,r) 	fp_copy(r,a)
# define fp_is_one(a)		bn_is_one(a,8)
# define fp_mul(r,a,b)		bn_barrett_mod_mul(r,a,b,BN254_P,BN254_P_MU,8)
# define fp_sqr(r,a)		bn_barrett_mod_sqr(r,a,BN254_P,BN254_P_MU,8)
#endif

void fp_tri(fp_t r, const fp_t a);
void fp_exp(fp_t r, const fp_t a, const uint32_t *e, int elen);
void fp_inv(fp_t r, const fp_t a);


typedef struct {
	fp_t X;
	fp_t Y;
	fp_t Z;
	int is_at_infinity;
} point_t;

typedef struct {
	fp_t X;
	fp_t Y;
} affine_point_t;

extern const point_t BN254_G1;


void point_set_infinity(point_t *R);
int point_is_at_infinity(const point_t *P);
void point_copy(point_t *R, const point_t *P);
void point_set_xy(point_t *R, const uint32_t x[8], const uint32_t y[8]);
void point_get_xy(const point_t *P, uint32_t x[8], uint32_t y[8]);
void point_get_affine(const point_t *P, affine_point_t *R);
void point_from_hex(point_t *R, const char *hex);
int point_equ_hex(const point_t *P, const char *hex);
void point_print(const char *str, const point_t *P);
void affine_point_print(const char *str, const affine_point_t *P);
void point_to_bytes(const point_t *P, uint8_t out[64]);
void point_from_bytes(point_t *P, const uint8_t in[64]);
void point_dbl(point_t *R, const point_t *P);
void point_add_affine(point_t *R, const point_t *P, const affine_point_t *Q);
void point_add_jacobian(point_t *R, const point_t *P, const point_t *Q);
void point_neg(point_t *R, const point_t *P);
#define point_add(R,P,Q) point_add_jacobian(R,P,Q)
void point_sub(point_t *R, const point_t *P, const point_t *Q);
void point_mul_affine_non_const_time(point_t *R, const uint32_t a[8], const affine_point_t *P);
void point_mul_affine_const_time(point_t *R, const uint32_t a[8], const affine_point_t *P);
#ifdef BN254_CONST_TIME
# define point_mul_affine(R,a,P)	point_mul_affine_const_time(R,a,P)
#else
# define point_mul_affine(R,a,P)	point_mul_affine_non_const_time(R,a,P)
#endif
void point_mul_generator(point_t *R, const uint32_t a[8]);
void point_multi_mul(point_t *R, const uint32_t a[8][8], const affine_point_t *P);
void point_multi_mul_affine_pre_compute(const point_t *P, point_t *T);
void point_multi_mul_affine_with_pre_compute(point_t *R, const uint32_t a[8][8], const point_t *T);
void point_multi_exponent(point_t *R, const fn_t *a, int alen, const affine_point_t *CRS);
void plonk_generate_crs(affine_point_t *CRS, const fn_t s, int n);


typedef fn_t *poly_t;

poly_t poly_new(int len);
void poly_free(poly_t *a);
void poly_copy(poly_t r, int *rlen, const poly_t a, int alen);
void poly_rand(poly_t r, int n);
int  poly_to_file(const poly_t a, int n, FILE *fp);
int  poly_from_file(poly_t r, int n, FILE *fp);
void poly_print(const char *str, const poly_t a, int alen);
void poly_add(poly_t r, int *rlen, const poly_t a, int alen, const poly_t b, int blen);
void poly_sub(poly_t r, int *rlen, const poly_t a, int alen, const poly_t b, int blen);
void poly_mul(poly_t r, int *rlen, const poly_t a, int alen, const poly_t b, int blen);
void poly_add_scalar(poly_t r, int *rlen, const poly_t a, int alen, const fn_t scalar);
void poly_sub_scalar(poly_t r, int *rlen, const poly_t a, int alen, const fn_t scalar);
void poly_mul_scalar(poly_t r, int *rlen, const poly_t a, int alen, const fn_t scalar);
void poly_mul_vector(poly_t r, int *rlen, const poly_t a, int alen, const fn_t *vec, int veclen);
void poly_add_blind(poly_t a, int *rlen, const poly_t b, int blen);
void poly_div_x_sub_scalar(poly_t a, int *alen, const fn_t scalar);
void poly_div_ZH(poly_t r, int *rlen, const poly_t a, int alen, int n);
void poly_eval(fn_t r, const fn_t *a, int alen, const fn_t x);
void poly_eval_L1(fn_t r, int n, const fn_t x);

int reverse_bits(int i, int nbits);
void fft(fn_t *vals, const fn_t *H, int nbits);
void poly_interpolate(fn_t *vals, const fn_t *H, int nbits);

#ifdef __cplusplus
}
#endif
#endif
