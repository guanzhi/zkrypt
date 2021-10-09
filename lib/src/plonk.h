#ifndef ZCRYPT_PLONK_H
#define ZCRYPT_PLONK_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif


// for testing only
void plonk_generate_crs(affine_point_t *CRS, const fn_t s, int n);

void plonk_generate_domain(fn_t *H, const fn_t w, int n);

int plonk_proof(
	const affine_point_t *CRS,
	const fn_t *H,
	const int nbits,
	const fn_t k1,
	const fn_t k2,
	const poly_t L1,
	const poly_t qM,
	const poly_t qL,
	const poly_t qR,
	const poly_t qO,
	const poly_t qC,
	const poly_t PI,
	const int *sigma,
	const poly_t S_sigma1,
	const poly_t S_sigma2,
	const poly_t S_sigma3,
	const fn_t *in,
	point_t *A,
	point_t *B,
	point_t *C,
	point_t *Z,
	point_t *T_lo,
	point_t *T_mid,
	point_t *T_hi,
	point_t *W_xi,
	point_t *W_xi_w,
	fn_t _a,
	fn_t _b,
	fn_t _c,
	fn_t _s_sigma1,
	fn_t _s_sigma2,
	fn_t _t,
	fn_t _r,
	fn_t _z_w);

#ifdef __cplusplus
}
#endif
#endif
