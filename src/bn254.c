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
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "bn254.h"


/*
BN254 curve parameters
from "Faster Explicit Formulas for Computing Pairings over Ordinary Curves"

u = -1 * (2^62 + 2^55 + 1)
p = 36 * u^4 + 36 * u^3 + 24 * u^2 + 6 * u + 1
n = 36 * u^4 + 36 * u^3 + 18 * u^2 + 6 * u + 1
t = 6 * u^2 + 1
*/


const uint32_t BN254_ONE[] = { 1,0,0,0,0,0,0,0 };

const uint32_t BN254_TWO[] = { 2,0,0,0,0,0,0,0 };


void bn_print(const char *str, const uint32_t *a, int k)
{
	printf("%s", str);
	while (k-- > 0) {
		printf("%08x", a[k]);
	}
	printf("\n");
}


// re order a and k, r should be together with k
void bn_set_word(uint32_t *r, uint32_t a, int k)
{
	r[0] = a;
	while (k-- > 1) {
		r[k] = 0;
	}
}

void bn_copy(uint32_t *r, const uint32_t *a, int k)
{
	while (k-- > 0) {
		r[k] = a[k];
	}
}

int bn_cmp(const uint32_t *a, const uint32_t *b, int k)
{
	while (k-- > 0) {
		if (a[k] > b[k]) return 1;
		else if (a[k] < b[k]) return -1;
	}
	return 0;
}

int bn_is_zero(const uint32_t *a, int k)
{
	while (k-- > 0) {
		if (a[k]) {
			return 0;
		}
	}
	return 1;
}

int bn_is_one(const uint32_t *a, int k)
{
	if (a[0] != 1) {
		return 0;
	}
	while (k-- > 1) {
		if (a[k]) {
			return 0;
		}
	}
	return 1;
}

int bn_add(uint32_t *r, const uint32_t *a, const uint32_t *b, int k)
{
	uint64_t w = 0;
	int i;
	for (i = 0; i < k; i++) {
		w += (uint64_t)a[i] + (uint64_t)b[i];
		r[i] = w & 0xffffffff;
		w >>= 32;
	}
	return (int)w;
}

int bn_sub(uint32_t *r, const uint32_t *a, const uint32_t *b, int k)
{
	int64_t w = 0;
	int i;
	for (i = 0; i < k; i++) {
		w += (int64_t)a[i] - (int64_t)b[i];
		r[i] = w & 0xffffffff;
		w >>= 32;
	}
	return (int)w;
}

void bn_mul(uint32_t *r, const uint32_t *a, const uint32_t *b, int k)
{
	uint64_t w;
	int i, j;
	for (i = 0; i < k; i++) {
		r[i] = 0;
	}
	for (i = 0; i < k; i++) {
		w = 0;
		for (j = 0; j < k; j++) {
			w += (uint64_t)r[i + j] + (uint64_t)a[i] * (uint64_t)b[j];
			r[i + j] = w & 0xffffffff;
			w >>= 32;
		}
		r[i + k] = w;
	}
}

void bn_mul_lo(uint32_t *r, const uint32_t *a, const uint32_t *b, int k)
{
	uint64_t w;
	int i, j;
	for (i = 0; i < k; i++) {
		r[i] = 0;
	}
	for (i = 0; i < k; i++) {
		w = 0;
		for (j = 0; j < k - i; j++) {
			w += (uint64_t)r[i + j] + (uint64_t)a[i] * (uint64_t)b[j];
			r[i + j] = w & 0xffffffff;
			w >>= 32;
		}
	}
}

static void uint32_to_bytes(uint32_t a, uint8_t out[4])
{
	out[0] = (a >> 24) & 0xff;
	out[1] = (a >> 16) & 0xff;
	out[2] = (a >> 8) & 0xff;
	out[3] = a & 0xff;
}

static uint32_t uint32_from_bytes(const uint8_t in[4])
{
	return ((uint32_t)in[0] << 24)
		| ((uint32_t)in[1] << 16)
		| ((uint32_t)in[2] << 8)
		| in[3];
}

void bn_to_bytes(const uint32_t *a, int k, uint8_t *out)
{
	while (k-- > 0) {
		uint32_to_bytes(a[k], out);
		out += 4;
	}
}

void bn_from_bytes(uint32_t *a, int k, const uint8_t *in)
{
	while (k-- > 0) {
		a[k] = uint32_from_bytes(in);
		in += 4;
	}
}

void bn_to_hex(const uint32_t *a, int k, char *out)
{
	while (k-- > 0) {
		sprintf(out, "%08x", a[k]);
		out += 8;
	}
	*out = 0;
}

static int hexchar2int(char c)
{
	if      ('0' <= c && c <= '9') return c - '0';
	else if ('a' <= c && c <= 'f') return c - 'a' + 10;
	else if ('A' <= c && c <= 'F') return c - 'A' + 10;
	else return -1;
}

static int hex2bin(const char *in, size_t inlen, uint8_t *out)
{
	int c;
	if (inlen % 2) {
		return -1;
	}

	while (inlen) {
		if ((c = hexchar2int(*in++)) < 0) {
			return -1;
		}
		*out = (uint8_t)c << 4;
		if ((c = hexchar2int(*in++)) < 0) {
			return -1;
		}
		*out |= (uint8_t)c;
		inlen -= 2;
		out++;
	}
	return 1;
}

int bn_from_hex(uint32_t *r, int k, const char *hex)
{
	uint8_t buf[k * 4];
	//assert(strlen(hex) == k * 8);
	if (hex2bin(hex, k * 8, buf) < 0) {
		return -1;
	}
	bn_from_bytes(r, k, buf);
	return 1;
}

int bn_equ_hex(const uint32_t *a, int k, const char *hex)
{
	uint32_t a_[k];
	bn_from_hex(a_, k, hex);
	return bn_cmp(a, a_, k) == 0;
}

int bn_rand(uint32_t *r, int k)
{
	FILE *fp = fopen("/dev/urandom", "rb");
	fread(r, 1, sizeof(uint32_t) * k, fp);
	fclose(fp);
	return 1;
}

// random r in [0, range-1]
int bn_rand_range(uint32_t *r, const uint32_t *range, int k)
{
	do {
		bn_rand(r, k);
	} while (bn_cmp(r, range, k) >= 0);
	return 1;
}




void bn_mod_add_non_const_time(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p, int k)
{
	bn_add(r, a, b, k);
	if (bn_cmp(r, p, k) >= 0) {
		bn_sub(r, r, p, k);
	}
}

void bn_mod_sub_non_const_time(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p, int k)
{
	if (bn_cmp(a, b, k) >= 0) {
		bn_sub(r, a, b, k);
	} else {
		bn_sub(r, b, a, k);
		bn_sub(r, p, r, k);
	}
}

void bn_mod_add_const_time(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p, int k)
{
	uint32_t t[2 * k];
	uint32_t *ret;
	int w;
	int i;

	bn_add(t, a, b, k);
	w = bn_sub(t + k, t, p, k);

	ret = t + k + k * w;
	for (i = 0; i < k; i++) {
		r[i] = ret[i];
	}
}

void bn_mod_sub_const_time(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p, int k)
{
	uint32_t r2[2 * k];
	uint32_t *rp;
	int w, i;

	w = bn_sub(r2 + k, a, b, k);
	bn_sub(r2, b, a, k);
	bn_sub(r2, p, r2, k);

	rp = r2 + k +  k * w;
	for (i = 0; i < k; i++) {
		r[i] = rp[i];
	}
}



void bn_mod_neg(uint32_t *r, const uint32_t *a, const uint32_t *p, int k)
{
	bn_sub(r, p, a, k);
}



static void bn_reduce_once(uint32_t *r, const uint32_t *a, const uint32_t *b, int k)
{
	uint32_t r2[k * 2]; // 为什么要这么大
	uint32_t *rp;
	int i;
	int w;

	// r0 = a
	for (i = 0; i < k; i++) {
		r2[i] = a[i];
	}

	// r1 = a - b
	w = bn_sub(r2 + k, a, b, k);

	rp = r2 + k + k * w;
	for (i = 0; i < k; i++) {
		r[i] = rp[i];
	}
}

void bn_barrett_mod_mul(uint32_t *r, // r[0..k] = a * b mod p
	const uint32_t *a, // a[0..k]
	const uint32_t *b, // b[0..k]
	const uint32_t *p, // p[0..k]
	const uint32_t *u, // u[0..k+1] = (2^32)^(2*k) // p, i.e. u = 2^512 // p
	int k)
{
	uint32_t z[2 * k];
	uint32_t q[2 * (k + 1)];
	uint32_t p_[k + 1];
	uint32_t t_[k + 1];
	uint32_t r_[k + 1];
	int i;

	for (i = 0; i < k; i++) {
		p_[i] = p[i];
	}
	p_[k] = 0;

	bn_mul(z, a, b, k);			// bn_print("z = a * b = ", z, 2 * k);
						// bn_print("zh = ", z + k - 1, k + 1);
						// bn_print("zl = ", z, k + 1);
	bn_mul(q, z + k - 1, u,	k + 1);		// bn_print("^q = ", q + k + 1, k + 1);
	bn_mul_lo(t_, q + k + 1, p_, k + 1);	// bn_print("q*n = ", t_, k + 1);
	bn_sub(r_, z, t_, k + 1);		// bn_print("r = ", r_, k + 1);

	bn_reduce_once(r_, r_, p_, k + 1);
	bn_reduce_once(r_, r_, p_, k);

	for (i = 0; i < k; i++) {
		r[i] = r_[i];
	}
}

void bn_barrett_mod_sqr(uint32_t *r, const uint32_t *a, const uint32_t *p,
	const uint32_t *u, int k)
{
	bn_barrett_mod_mul(r, a, a, p, u, k);
}

void bn_barrett_mod_exp(uint32_t *r, const uint32_t *a, const uint32_t *e, const uint32_t *p,
	const uint32_t *u, // u = mu(p) = 2^512 // p, len(u) = k + 1
	int k)
{
	uint32_t t[k];
	uint32_t w;
	int i, j;

	// t = 1
	bn_set_one(t, k);

	for (i = k - 1; i >= 0; i--) {
		w = e[i];
		for (j = 0; j < 32; j++) {
			bn_barrett_mod_sqr(t, t, p, u, k);
			// cuda support this branch
			if (w & 0x80000000) {
				bn_barrett_mod_mul(t, t, a, p, u, k);
			}
			w <<= 1;
		}
	}

	bn_copy(r, t, k);
}

void bn_barrett_mod_inv(uint32_t *r, const uint32_t *a, const uint32_t *p, const uint32_t *u, int k)
{
	uint32_t e[k];
	int i;

	e[0] = 2;
	for (i = 1; i < k; i++) {
		e[i] = 0;
	}

	bn_sub(e, p, e, k);
	bn_barrett_mod_exp(r, a, e, p, u, k);
}



// mont(aR, bR) = aR * bR * R^-1 = abR (mod p)
void bn_mont_mod_mul(uint32_t *r, const uint32_t *a, const uint32_t *b, const uint32_t *p,
	const uint32_t *p_inv_neg, int k)
{
	uint32_t z[k * 2];
	uint32_t c[k * 2];
	uint32_t t[k];
	int i;

	bn_mul(z, a, b, k);
	bn_mul_lo(t, z, p_inv_neg, k);
	bn_mul(c, t, p, k);
	bn_add(c, c, z, k * 2);
	if (bn_cmp(c + k, p, k) >= 0) {
		bn_sub(c + k, c + k, p, k);
	}

	for (i = 0; i < k; i++) {
		r[i] = c[k + i];
	}
}

void bn_mont_mod_sqr(uint32_t *r, const uint32_t *a, const uint32_t *p,
	const uint32_t *p_inv_neg, int k)
{
	bn_mont_mod_mul(r, a, a, p, p_inv_neg, k);
}

void bn_mont_mod_exp(
	uint32_t *r,
	const uint32_t *a,
	const uint32_t *e,
	const uint32_t *p,
	const uint32_t *p_inv_neg,
	int k)
{
	uint32_t t[k];
	uint32_t w;
	int i, j;

	bn_set_one(t, k);

	for (i = k - 1; i >= 0; i--) {
		w = e[i];
		for (j = 0; j < 32; j++) {
			bn_mont_mod_sqr(t, t, p, p_inv_neg, k);
			if (w & 0x80000000) {
				bn_mont_mod_mul(t, t, a, p, p_inv_neg, k);
			}
			w <<= 1;
		}
	}

	bn_copy(r, t, k);
}

void bn_mont_mod_inv(uint32_t *r, const uint32_t *a, const uint32_t *p,
	const uint32_t *p_inv_neg, int k)
{
	uint32_t e[k];
	int i;

	e[0] = 2;
	for (i = 1; i < k; i++) {
		e[i] = 0;
	}

	bn_sub(e, p, e, k);
	bn_mont_mod_exp(r, a, e, p, p_inv_neg, k);
}

// mont(a, R^2) = a * R^2 * R^-1 = a * R mod p
void bn_mont_set(uint32_t *r,
	const uint32_t *a,
	const uint32_t *R_sqr,
	const uint32_t *p,
	const uint32_t *p_inv_neg,
	int k)
{
	bn_mont_mod_mul(r, a, R_sqr, p, p_inv_neg, k);
}

// mont(aR, 1) = aR * 1 * R^-1 = a (mod p)
void bn_mont_get(uint32_t *r,
	const uint32_t *a,
	const uint32_t *p,
	const uint32_t *p_inv_neg,
	int k)
{
	uint32_t one[k];
	bn_set_one(one, k);
	bn_mont_mod_mul(r, a, one, p, p_inv_neg, k);
}

// 这个函数仅用于测试
void bn_mod_mul_montgomery(uint32_t *r, const uint32_t *a, const uint32_t *b,
	const uint32_t *R_sqr, const uint32_t *p, const uint32_t *p_inv_neg, int k)
{
	uint32_t a_[k];
	uint32_t b_[k];
	uint32_t r_[k];

	bn_mont_set(a_, a, R_sqr, p, p_inv_neg, k);
	bn_mont_set(b_, b, R_sqr, p, p_inv_neg, k);
	bn_mont_mod_mul(r_, a_, b_, p, p_inv_neg, k);
	bn_mont_get(r, r_, p, p_inv_neg, k);
}

void bn_mod_sqr_montgomery(uint32_t *r, const uint32_t *a,
	const uint32_t *R_sqr, const uint32_t *p, const uint32_t *p_inv_neg, int k)
{
	uint32_t a_[k];
	uint32_t r_[k];

	bn_mont_set(a_, a, R_sqr, p, p_inv_neg, k);
	bn_mont_mod_sqr(r_, a_, p, p_inv_neg, k);
	bn_mont_get(r, r_, p, p_inv_neg, k);
}








void fn_tri(fn_t r, const fn_t a)
{
	fn_t t;
	fn_dbl(t, a);
	fn_add(r, t, a);
}

void fn_exp(fn_t r, const fn_t a, const uint32_t *e, int elen)
{
	uint32_t t[8];
	uint32_t w;
	int i, j;

	fn_set_one(t);

	for (i = elen - 1; i >= 0; i--) {
		w = e[i];
		for (j = 0; j < 32; j++) {
			fn_sqr(t, t);
			if (w & 0x80000000) {
				fn_mul(t, t, a);
			}
			w <<= 1;
		}
	}

	fn_copy(r, t);
}





void fn_inv(uint32_t *r, const uint32_t *a)
{
	uint32_t e[8];
	int i;

	e[0] = 2;
	for (i = 1; i < 8; i++) {
		e[i] = 0;
	}

	bn_sub(e, BN254_N, e, 8);
	fn_exp(r, a, e, 8);
}





void fp_tri(fp_t r, const fp_t a)
{
	fp_t t;
	fp_dbl(t, a);
	fp_add(r, t, a);
}

void fp_exp(fp_t r, const fp_t a, const uint32_t *e, int elen)
{
	uint32_t t[8];
	uint32_t w;
	int i, j;

	fp_set_one(t);

	for (i = elen - 1; i >= 0; i--) {
		w = e[i];
		for (j = 0; j < 32; j++) {
			fp_sqr(t, t);
			if (w & 0x80000000) {
				fp_mul(t, t, a);
			}
			w <<= 1;
		}
	}

	fp_copy(r, t);
}

#if 1
void fp_inv(uint32_t *r, const uint32_t *a)
{
	uint32_t e[8];
	int i;

	e[0] = 2;
	for (i = 1; i < 8; i++) {
		e[i] = 0;
	}

	bn_sub(e, BN254_P, e, 8);
	fp_exp(r, a, e, 8);
}
#else
// 这个实际上是和p相关的，我们应该检查一下p和默认的p是不是一样的
// a^-1 = a^(p-2) (mod p)
void fp_inv(fp_t r, const fp_t a)
{
	int i;
	fp_copy(r, a);
	fp_sqr(r, r);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 3; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 5; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 29; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	for (i = 0; i < 3; i++) {
		fp_sqr(r, r); fp_mul(r, r, a);
	}
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 3; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 3; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 35; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 4; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 4; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 4; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 43; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r);
	for (i = 0; i < 3; i++) {
		fp_sqr(r, r); fp_mul(r, r, a);
	}
	fp_sqr(r, r);
	fp_sqr(r, r); fp_mul(r, r, a);
	fp_sqr(r, r);
	fp_sqr(r, r);
	for (i = 0; i < 3; i++) {
		fp_sqr(r, r); fp_mul(r, r, a);
	}
	for (i = 0; i < 51; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
	for (i = 0; i < 3; i++) {
		fp_sqr(r, r);
	}
	fp_sqr(r, r); fp_mul(r, r, a);
}
#endif

void point_set_infinity(point_t *R)
{
	fp_set_one(R->X);
	fp_set_one(R->Y);
	fp_set_zero(R->Z);
	R->is_at_infinity = 1;
}

int point_is_at_infinity(const point_t *P)
{
	return P->is_at_infinity;
}

// P is affine ?
void point_copy(point_t *R, const point_t *P)
{
	fp_copy(R->X, P->X);
	fp_copy(R->Y, P->Y);
	fp_copy(R->Z, P->Z);
	R->is_at_infinity = P->is_at_infinity;
}

// x, y should not be in Montgomery field
void point_set_xy(point_t *R, const uint32_t x[8], const uint32_t y[8])
{
	fp_set_bn(R->X, x);
	fp_set_bn(R->Y, y);
	fp_set_one(R->Z);
}

#if 0
void point_get_xy(const point_t *P, uint32_t x[8], uint32_t y[8])
{
	fp_t z_inv;

	if (fp_is_one(P->Z)) {
		fp_get_bn(P->X, x);
		fp_get_bn(P->Y, y);
	} else {
		fp_inv(z_inv, P->Z);
		if (y) {
			fp_mul(y, P->Y, z_inv);
		}
		fp_sqr(z_inv, z_inv);
		fp_mul(x, P->X, z_inv);
		fp_get_bn(x, x);
		if (y) {
			fp_mul(y, y, z_inv);
			fp_get_bn(y, y);
		}
	}
}
#else
void point_get_xy(const point_t *P, uint32_t x[8], uint32_t y[8])
{
	fp_t Z_inv;

	fp_inv(Z_inv, P->Z);
	fp_mul(y, P->Y, Z_inv);
	fp_sqr(Z_inv, Z_inv);
	fp_mul(x, P->X, Z_inv);
	fp_get_bn(x, x);
	fp_mul(y, y, Z_inv);
	fp_get_bn(y, y);
}
#endif

void point_get_affine(const point_t *P, affine_point_t *R)
{
	fp_t Z_inv;

	fp_inv(Z_inv, P->Z);
	fp_mul(R->Y, P->Y, Z_inv);
	fp_sqr(Z_inv, Z_inv);
	fp_mul(R->X, P->X, Z_inv);
	//fp_get_bn(R->X, R->X);
	fp_mul(R->Y, R->Y, Z_inv);
	//fp_get_bn(R->Y, R->Y);
}

void point_from_hex(point_t *R, const char *hex)
{
	fp_from_hex(R->X, hex);
	fp_from_hex(R->Y, hex + 64);
	fp_set_one(R->Z);
	R->is_at_infinity = 0;
}

int point_equ_hex(const point_t *P, const char *hex)
{
	return 1;
}


void point_print(const char *str, const point_t *P)
{
	uint32_t x[8];
	uint32_t y[8];
	int i;

	point_get_xy(P, x, y);

	if (P->is_at_infinity) {
		printf("infinity\n");
	}

	printf("%s", str);
	for (i = 7; i >=0; i--) {
		printf("%08x", x[i]);
	}
	printf(" ");
	for (i = 7; i >= 0; i--) {
		printf("%08x", y[i]);
	}
	printf("\n");
}

void affine_point_print(const char *str, const affine_point_t *P)
{
	int i;
	printf("%s", str);
	for (i = 7; i >= 0; i--) {
		printf("%08x", P->X[i]);
	}
	printf(" ");
	for (i = 7; i >=0; i--) {
		printf("%08x", P->Y[i]);
	}
	printf("\n");
}

void point_to_bytes(const point_t *P, uint8_t out[64])
{
	uint32_t x[8];
	uint32_t y[8];
	point_get_xy(P, x, y);
	bn_to_bytes(x, 8, out);
	bn_to_bytes(y, 8, out + 32);
}

void point_from_bytes(point_t *P, const uint8_t in[64])
{
	fp_from_bytes(P->X, in);
	fp_from_bytes(P->Y, in + 32);
	fp_set_one(P->Z);
}

// 支持 P = (1:1:0)
// P 必须为Jacobian格式
void point_dbl(point_t *R, const point_t *P)
{
	fp_t T_0;
	fp_t T_1;
	fp_t T_2;
	fp_t T_3;
	fp_t T_4;

	fp_sqr(T_0, P->X);
	fp_tri(T_0, T_0);
	fp_sqr(T_1, T_0);
	fp_sqr(T_2, P->Y);
	fp_mul(T_3, P->X, T_2);
	fp_dbl(T_3, T_3);
	fp_dbl(T_3, T_3);
	fp_dbl(T_4, T_3);
	fp_sub(T_1, T_1, T_4);
	fp_sub(T_3, T_3, T_1);
	fp_mul(T_0, T_0, T_3);
	fp_dbl(T_2, T_2);
	fp_sqr(T_2, T_2);
	fp_dbl(T_2, T_2);
	fp_sub(T_0, T_0, T_2);
	fp_mul(T_2, P->Y, P->Z);
	fp_dbl(T_2, T_2);

	fp_copy(R->X, T_1);
	fp_copy(R->Y, T_0);
	fp_copy(R->Z, T_2);
	R->is_at_infinity = P->is_at_infinity;
}

//
// P 不可以是无穷远点，当 P = (1:1:0) 时，R只能得到(0:0:0)，因此调用方必须显式处理P==0的情况
// P 必须是Jacobian坐标
// Q 必须是仿射坐标，因此Q一定不是无穷远点

// 当 (X_1/Z_1^2, Y_1/Z_1^3) == (x_2, -y_2), 即 P + Q == 0 时，结果为(X3:Y3:0), 其中X3,Y3可能不为1，因此这个函数需要对其处理
// 点加计算过程是有可能产生一个O点的，一旦结果为0点，那么后续的计算就又成为特殊情况了，因此必须支持0点的计算！

// 当 P==Q 时，这个函数的结果是错的，因此调用方要保证 P != Q
void point_add_affine(point_t *R, const point_t *P, const affine_point_t *Q)
{
	fp_t T_0;
	fp_t T_1;
	fp_t T_2;
	fp_t T_3;
	fp_t T_4;
	fp_t T_5;

	if (P->is_at_infinity) {
		//point_copy(R, Q); // Q 是affine的时候正确吗？
		fp_copy(R->X, Q->X);
		fp_copy(R->Y, Q->Y);
		fp_set_one(R->Z);
		R->is_at_infinity = 0;
		return;
	}

	fp_sqr(T_0, P->Z);
	fp_mul(T_1, T_0, P->Z);
	fp_mul(T_1, Q->Y, T_1);
	fp_sub(T_1, T_1, P->Y);
	fp_mul(T_0, Q->X, T_0);
	fp_sub(T_2, T_0, P->X);
	fp_add(T_0, T_0, P->X);

	fp_sqr(T_3, T_2);
	fp_mul(T_4, T_3, T_0);
	fp_sqr(T_5, T_1);
	fp_sub(T_5, T_5, T_4);
	fp_mul(T_4, P->X, T_3);
	fp_sub(T_4, T_4, T_5);
	fp_mul(T_4, T_1, T_4);
	fp_mul(T_3, T_3, T_2);
	fp_mul(T_3, P->Y, T_3);
	fp_sub(T_4, T_4, T_3);
	fp_mul(T_3, T_2, P->Z);

	fp_copy(R->X, T_5);
	fp_copy(R->Y, T_4);
	fp_copy(R->Z, T_3);
	R->is_at_infinity = 0;

	if (fp_is_zero(T_3)) {
		point_set_infinity(R);
	}
}

// 可能出现 abG + baG 的情况，即需要调用倍点
// FIMXE 对特殊情况的处理
void point_add_jacobian(point_t *R, const point_t *P, const point_t *Q)
{
	fp_t T_1;
	fp_t T_2;
	fp_t T_3;
	fp_t T_4;
	fp_t T_5;
	fp_t T_6;
	fp_t T_7;
	fp_t T_8;

	if (P->is_at_infinity) {
		point_copy(R, Q);
		return;
	}

	fp_sqr(T_1, P->Z);	// T_1 = Z_1^2
	fp_sqr(T_2, Q->Z);	// T_2 = Z_2^2
	fp_mul(T_3, Q->X, T_1);	// T_3 = X_2 * Z_1^2
	fp_mul(T_4, P->X, T_2); // T_4 = X_1 * Z_2^2
	fp_add(T_5, T_3, T_4);	// T_5 = X_2 * Z_1^2 + X_1 * Z_2^2 = C
	fp_sub(T_3, T_3, T_4);	// T_3 = X_2 * Z_1^2 - X_1 * Z_2^2 = B
	fp_mul(T_1, T_1, P->Z);	// T_1 = Z_1^3
	fp_mul(T_1, T_1, Q->Y);	// T_1 = Y_2 * Z_1^3
	fp_mul(T_2, T_2, Q->Z);	// T_2 = Z_2^3
	fp_mul(T_2, T_2, P->Y);	// T_2 = Y_1 * Z_2^3
	fp_add(T_6, T_1, T_2);	// T_6 = Y_2 * Z_1^3 + Y_1 * Z_2^3 = D
	fp_sub(T_1, T_1, T_2);	// T_1 = Y_2 * Z_1^3 - Y_1 * Z_2^3 = A

	if (fp_is_zero(T_1) && fp_is_zero(T_3)) {
		point_dbl(R, P);
		return;
	}

	if (fp_is_one(T_1) && fp_is_zero(T_6)) {
		point_set_infinity(R);
		return;
	}

	fp_sqr(T_6, T_1);	// T_6 = A^2
	fp_mul(T_7, T_3, P->Z);	// T_7 = B * Z_1
	fp_mul(T_7, T_7, Q->Z);	// T_7 = B * Z_1 * Z_2 = Z_3
	fp_sqr(T_8, T_3);	// T_8 = B^2
	fp_mul(T_5, T_5, T_8);	// T_5 = B^2 * C
	fp_mul(T_3, T_3, T_8);	// T_3 = B^3
	fp_mul(T_4, T_4, T_8);	// T_4 = B^2 * X_1 * Z_2^2
	fp_sub(T_6, T_6, T_5);	// T_6 = A^2 - B^2 * C = X_3
	fp_sub(T_4, T_4, T_6);	// T_4 = B^2 * X_1 * Z_2^2 - X_3
	fp_mul(T_1, T_1, T_4);	// T_1 = A * (B^2 * X_1 * Z_2^2 - X_3)
	fp_mul(T_2, T_2, T_3);	// T_2 = B^3 * Y_1 * Z_1^3
	fp_sub(T_1, T_1, T_2);	// T_1 = A * (B^2 * X_1 * Z_2^2 - X_3) - B^3 * Y_1 * Z_1^3 = Y_3

	fp_copy(R->X, T_6);
	fp_copy(R->Y, T_1);
	fp_copy(R->Z, T_7);
	R->is_at_infinity = 0;
}

#define point_add(R,P,Q) point_add_jacobian(R,P,Q)

void point_neg(point_t *R, const point_t *P)
{
	fp_copy(R->X, P->X);
	fp_neg(R->Y, P->Y);
	fp_copy(R->Z, P->Z);
}

// 这里面Q是Affine/Jacobina??
void point_sub(point_t *R, const point_t *P, const point_t *Q)
{
	point_t _T, *T = &_T;
	point_neg(T, Q);
	point_add(R, P, T);
}

// point_add要求输入的P必须为仿射坐标
// 这个函数应该改写一下，保证a的值是从1开始的
// 如果P != O, 那么我们应该从a == 1开始这个循环
// 有没有必要将P设置位Affine Coordinates?

// 实际上点乘函数是一个特殊情况，在点乘中，每次point_add的输入Q都是一个Affine坐标的点
void point_mul_affine_non_const_time(point_t *R, const uint32_t a[8], const affine_point_t *P)
{
	uint32_t bits;
	int nbits;
	int k;

	point_set_infinity(R);

	bits = a[7] << 2;
	nbits = 30;
	while (nbits-- > 0) {
		point_dbl(R, R);
		if (bits & 0x80000000) {
			point_add_affine(R, R, P);
		}
		bits <<= 1;
	}

	for (k = 6; k >= 0; k--) {
		bits = a[k];
		nbits = 32;
		while (nbits-- > 0) {
			point_dbl(R, R);
			if (bits & 0x80000000) {
				point_add_affine(R, R, P);
			}
			bits <<= 1;
		}
	}
}

// 标量乘法和指数运算时，这个运算数是一个标量，因此必须是一个整数，否则Mont格式可能产生错误
void point_mul_affine_const_time(point_t *R, const uint32_t a[8], const affine_point_t *P)
{
	uint32_t bits;
	int nbits;
	int k;


	point_set_infinity(R);

	bits = a[7] << 2;
	nbits = 30;
	while (nbits-- > 0) {
		point_t Rs[2];
		int do_add = (bits & 0x80000000) >> 31;
		point_dbl(&Rs[0], R);
		point_add_affine(&Rs[1], &Rs[0], P);
		point_copy(R, &Rs[do_add]);
		bits <<= 1;
	}

	for (k = 6; k >= 0; k--) {
		bits = a[k];
		nbits = 32;
		while (nbits-- > 0) {
			point_t Rs[2];
			int do_add = (bits & 0x80000000) >> 31;
			point_dbl(&Rs[0], R);
			point_add_affine(&Rs[1], &Rs[0], P);
			point_copy(R, &Rs[do_add]);
			bits <<= 1;
		}
	}
}



void point_mul_generator(point_t *R, const uint32_t a[8])
{
	point_mul_affine(R, a, (affine_point_t *)&BN254_G1);
}


// 有可能出现O点
// a1b1*G1 + a2b2*G1 = (a1b1 + a2b2)G1 ==0 => a1b1 + a2b2 = n
// multi_mul中，P[]占用非常多的存储空间，因此可以考虑引入 affine_point, jacobian_point 两种不同的类型
// affine_point 是一种特殊情况
void point_multi_mul(point_t *R, const uint32_t a[8][8], const affine_point_t *P)
{
	point_t _T, *T = &_T;
	int i;
	point_set_infinity(R);

	for (i = 0; i < 8; i++) {
		point_mul_affine(T, a[i], &P[i]); // 这里的P[i]是输入的固定点，总是Affine的
		point_add(R, R, T);  // 显然在这个场景中，T点是Jacobian Coordinates的,可能出现R == T的倍点情况
	}
}


// 当point坐标是Mont时，预计算的表T中的也必须是Mont值
// 应该把这个表整个存储到内存中

void point_multi_mul_affine_pre_compute(const point_t *P, point_t *T)
{
	int i;
	for (i = 0; i < 255; i++) {
		uint8_t iv = i + 1;
		int k;
		point_set_infinity(&T[i]);
		for (k = 0; k < 8; k++) {
			if (iv & 0x01) {
				point_add(&T[i], &T[i], &P[k]);
			}
			iv >>= 1;
		}
	}
}

// 这里一个问题，是否 2R == T[i]
// 这里 T[i] = a0 G + a1 w G + a2 w^2 G + ...

void point_multi_mul_affine_with_pre_compute(point_t *R, const uint32_t a[8][8], const point_t *T)
{
	uint32_t nwords;
	uint32_t bits[8];
	uint32_t nbits;
	uint8_t u;
	int i;

	point_set_infinity(R);

	for (i = 0; i < 8; i++) {
		bits[i] = (a[i][7]) << 2;
	}
	nbits = 30;
	while (nbits--) {
		u = 0;
		for (i = 7; i >= 0; i--) {
			u |= ((bits[i] & 0x80000000) >> 31) << i;
		}
		point_dbl(R, R);
		if (u) {
			point_add(R, R, &T[u-1]);
		}
		for (i = 0; i < 8; i++) {
			bits[i] <<= 1;
		}
	}

	nwords = 7;
	while (nwords--) {
		for (i = 0; i < 8; i++) {
			bits[i] = a[i][nwords];
		}
		nbits = 32;
		while (nbits--) {
			uint8_t u = 0;
			for (i = 7; i >= 0; i--) {
				u |= ((bits[i] & 0x80000000) >> 31) << i;
			}
			point_dbl(R, R);
			if (u) {
				point_add(R, R, &T[u-1]);
			}
			for (i = 0; i < 8; i++) {
				bits[i] <<= 1;
			}

		}
	}
}

void point_multi_exponent(point_t *R, const fn_t *a, int alen, const affine_point_t *CRS)
{
	point_t T;
	int i;

	point_set_infinity(R);
	for (i = 0; i < alen; i++) {
		point_mul_affine(&T, a[i], &CRS[i]);
		point_add(R, R, &T);
	}
}


poly_t poly_new(int len)
{
	fn_t *ret;
	if (!(ret = (fn_t *)malloc(sizeof(fn_t) * len))) {
		fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
		return NULL;
	}
	memset(ret, 0, sizeof(fn_t) * len);
	return ret;
}

void poly_free(poly_t *a)
{
	free(a);
}

void poly_copy(poly_t r, int *rlen, const poly_t a, int alen)
{
	int i;
	for (i = 0; i < alen; i++) {
		fn_copy(r[i], a[i]);
	}
	*rlen = alen;
}

void poly_rand(poly_t r, int n)
{
	int i;
	for (i = 0; i < n; i++) {
		fn_rand(r[i]);
	}
}

int poly_to_file(const poly_t a, int n, FILE *fp)
{
	char buf[128];
	int i;

	for (i = 0; i < n; i++) {
		fn_to_hex(a[i], buf);
		fprintf(fp, "%s\n", buf);
	}
	return 1;
}

int poly_from_file(poly_t r, int n, FILE *fp)
{
	char buf[128];
	int i;

	for (i = 0; i < n; i++) {
		if (!fgets(buf, sizeof(buf), fp)) {
			fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
			return -1;
		}
		if (strlen(buf) < 64) {
			fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
			return -1;
		}
		fn_from_hex(r[i], buf);
	}
	return 1;
}

void poly_print(const char *str, const poly_t a, int alen)
{
	int i;
	printf("%s\n", str);
	for (i = 0; i < alen; i++) {
		printf("  %d", i);
		bn_print(" ", a[i], 8);
	}
	printf("\n");
}

void poly_add(poly_t r, int *rlen, const poly_t a, int alen, const poly_t b, int blen)
{
	const fn_t *fa = a;
	const fn_t *fb = b;
	int falen = alen;
	int fblen = blen;
	int i;

	if (alen > blen) {
		fa = b;
		fb = a;
		falen = blen;
		fblen = alen;
	}

	for (i = 0; i < falen; i++) {
		fn_add(r[i], fa[i], fb[i]);
	}
	for (; i < fblen; i++) {
		fn_copy(r[i], fb[i]);
	}

	*rlen = fblen;
}

void poly_sub(poly_t r, int *rlen, const poly_t a, int alen, const poly_t b, int blen)
{
	int i;

	if (alen >= blen) {
		for (i = 0; i < blen; i++) {
			fn_sub(r[i], a[i], b[i]);
		}
		for (; i < alen; i++) {
			fn_copy(r[i], a[i]);
		}
		*rlen = alen;
	} else {
		for (i = 0; i < alen; i++) {
			fn_sub(r[i], a[i], b[i]);
		}
		for (; i < blen; i++) {
			fn_neg(r[i], b[i]);
		}
		*rlen = blen;
	}
}

void poly_mul(poly_t r, int *rlen, const poly_t a, int alen, const poly_t b, int blen) // TODO: poly_mul_ex(buffer); # foo foo_ex(buf);
{
	int i, j;
	fn_t t;
	fn_t *tmp = poly_new(alen + blen - 1);

	for (j = 0; j < blen; j++) {
		fn_mul(tmp[j], a[0], b[j]);
	}
	for (i = 1; i < alen; i++) {
		for (j = 0; j < blen - 1; j++) {
			fn_mul(t, a[i], b[j]);
			fn_add(tmp[i + j], tmp[i + j], t);
		}
		fn_mul(tmp[i + j], a[i], b[j]);
	}
	*rlen = alen + blen - 1;
	poly_copy(r, rlen, tmp, alen + blen - 1);
}

void poly_add_scalar(poly_t r, int *rlen, const poly_t a, int alen, const fn_t scalar)
{
	fn_add(r[0], a[0], scalar);
	if (r != a) {
		int i;
		for (i = 1; i < alen; i++) {
			fn_copy(r[i], a[i]);
		}
	}
	*rlen = alen;
}

void poly_sub_scalar(poly_t r, int *rlen, const poly_t a, int alen, const fn_t scalar)
{
	fn_sub(r[0], a[0], scalar);
	if (r != a) {
		int i;
		for (i = 1; i < alen; i++) {
			fn_copy(r[i], a[i]);
		}
	}
	*rlen = alen;
}

void poly_mul_scalar(poly_t r, int *rlen, const poly_t a, int alen, const fn_t scalar)
{
	int i;
	for (i = 0; i < alen; i++) {
		fn_mul(r[i], a[i], scalar);
	}
	*rlen = alen;
}

// r[i] = a[i] * vec[i]
void poly_mul_vector(poly_t r, int *rlen, const poly_t a, int alen, const fn_t *vec, int veclen)
{
	int i;
	for (i = 0; i < alen; i++) {
		fn_mul(r[i], a[i], vec[i % veclen]);
	}
	*rlen = alen;
}

// r(x) = a(x) + (b_0 + b_1 * x [+ b_2 * x^2]) * (x^n - 1), n = len(a)
void poly_add_blind(poly_t a, int *rlen, const poly_t b, int blen)
{
	int i;
	for (i = 0; i < blen; i++) {
		fn_sub(a[i], a[i], b[i]);
		fn_copy(a[*rlen + i], b[i]);
	}
	(*rlen) += blen;
}

// 未经完全测试
// a(x) = a(x)/(x - scalar)
void poly_div_x_sub_scalar(poly_t a, int *alen, const fn_t scalar)
{
	int i;
	fn_t t;
	for (i = *alen - 2; i >= 0; i--) {
		fn_mul(t, a[i + 1], scalar);
		fn_add(a[i], a[i], t);
	}
	// a[0] 是余数，并且应该为0
	// TODO: check if a[0] == 0
	for (i = 0; i < *alen - 1; i++) {
		fn_copy(a[i], a[i + 1]);
	}
	(*alen)--;
}

// 未经测试，注意a(x)要足够长
// r(x) = a(x)/(x^n - 1)
void poly_div_ZH(poly_t r, int *rlen, const poly_t a, int alen, int n)
{
	int i;
	for (i = 0; i < n; i++) {
		fn_neg(r[i], a[i]);
	}
	for (; i < alen; i++) {
		fn_sub(r[i], r[i - n], a[i]);
	}
	*rlen = alen - n;
}

void poly_eval(fn_t r, const fn_t *a, int alen, const fn_t x)
{
	fn_copy(r, a[--alen]);
	while (--alen >= 0) {
		fn_mul(r, r, x);
		fn_add(r, r, a[alen]);
	}
}

// 需要优化吗？
void fn_mul_word(fn_t r, const fn_t a, uint32_t word)
{
	fn_t b;
	bn_set_word(b, word, 8);
	fn_mul(r, a, b);
}

// L1(x) = (x^n - 1)/(n * (x - 1))
void poly_eval_L1(fn_t r, int n, const fn_t x)
{
	fn_t t;
	fn_exp(t, x, (uint32_t *)&n, 1);
	fn_sub(t, t, BN254_ONE);

	fn_sub(r, x, BN254_ONE);
	fn_mul_word(r, r, n);
	fn_inv(r, r);
	fn_mul(r, r, t);
}

int reverse_bits(int i, int nbits)
{
	unsigned int r = 0;
	while (nbits-- > 0) {
		r = (r << 1) | (i & 1);
		i >>= 1;
	}
	return r;
}

void fft(fn_t *vals, const fn_t *H, int nbits)
{
	int n = 1 << nbits;
	int Hoffset = n/2;
	int i, j, k;

	for (i = 0; i < n; i++) {
		fn_t tmp;
		int i_rev = reverse_bits(i, nbits);
		if (i < i_rev) {
			fn_copy(tmp, vals[i]);
			fn_copy(vals[i], vals[i_rev]);
			fn_copy(vals[i_rev], tmp);
		}
	}

	for (i = 0; i < nbits; i++) {

		int half_len = 1 << i;
		int full_len = half_len << 1;
		int count = (1 << nbits)/full_len;

		for (j = 0; j < count; j++) {
			int Li = full_len * j;
			int Ri = Li + half_len;
			int Hi = 0;

			for (k = 0; k < half_len; k++) {
				fn_t x;
				fn_t y;

				fn_copy(x, vals[Li]);
				fn_copy(y, vals[Ri]);
				fn_mul(y, y, H[Hi]);
				fn_add(vals[Li], x, y);
				fn_sub(vals[Ri], x, y);

				Li++;
				Ri++;
				Hi += Hoffset;
			}
		}
		Hoffset /= 2;
	}
}

void poly_interpolate(fn_t *vals, const fn_t *H, int nbits)
{
	int n;
	fn_t n_inv;
	int i;

	n = 1 << nbits;

	fft(vals, H, nbits);

	bn_set_word(n_inv, n, 8);



	fn_inv(n_inv, n_inv);


	poly_mul_scalar(vals, &n, vals, n, n_inv);

	for (i = 1; i < (n + 1)/2; i++) {
		fn_copy(n_inv, vals[i]);
		fn_copy(vals[i], vals[n - i]);
		fn_copy(vals[n - i], n_inv);
	}
}
