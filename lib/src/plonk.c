#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "sha2.h"
#include "bn254.h"
#include "plonk.h"
#include "bn254_test.h"

void plonk_generate_crs(affine_point_t *CRS, const fn_t s, int n)
{
	fn_t s_exp;
	point_t P;
	int i;

	fn_set_word(s_exp, 1);
	for (i = 0; i < n; i++) {
		point_mul_generator(&P, s_exp);
		point_get_affine(&P, &CRS[i]);
		fn_mul(s_exp, s_exp, s);
	}
}

// generate domain[0..n-1] = H, domain[n..2n-1] = k1*H, domain[2n..3n-1] = k2*H
void plonk_generate_domain(
	fn_t *domain,
	const fn_t w,
	//const fn_t k1,
	//const fn_t k2,
	int n)
{
	int i;
	fn_set_one(domain[0]);
	//fn_copy(domain[n], k1);
	//fn_copy(domain[n * 2], k2);
	for (i = 1; i < n; i++) {
		fn_mul(domain[i], domain[i - 1], w);
		//fn_mul(domain[n + i], domain[i], k1);
		//fn_mul(domain[n * 2 + i], domain[i], k2);
	}
}

int plonk_proof(
	const affine_point_t *CRS, // = [G1, [s]G1, [s^2]G1, ...]
	const fn_t *H, // = [1, w, w^2, ...]
	const int nbits, // len(qM) = len(qL) = ... = len(qC) = 2^nbits
	const fn_t k1,
	const fn_t k2,
	const poly_t L1,
	const poly_t qM,
	const poly_t qL,
	const poly_t qR,
	const poly_t qO,
	const poly_t qC,
	const poly_t PI,
	const int *sigma, // permutation of [0, 1, 2, ...], index of H, len(sigma) = 3 * n
	const poly_t S_sigma1,
	const poly_t S_sigma2,
	const poly_t S_sigma3,
	const fn_t *in, // len(in) = 3*n
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
	fn_t _z_w)
{

	int n = 1 << nbits;


	fn_t *a = poly_new(n + 2);
	fn_t *b = poly_new(n + 2);
	fn_t *c = poly_new(n + 2);
	fn_t *z = poly_new(n + 3);
	fn_t *t = poly_new(4 * n + 6);
	fn_t *r = poly_new(4 * n + 6);
	fn_t *g = poly_new(4 * n + 6);
	fn_t blind[9];

	fn_t fz[2]; // fz(X) = beta * X + gamma
	fn_t *f1 = poly_new(4 * n + 6);
	fn_t *f2 = poly_new(4 * n + 6);
	int f1_len = 0;
	int f2_len = 0;

	int a_len = 0;
	int b_len = 0;
	int c_len = 0;
	int z_len = 0;
	int t_len = 0;
	int r_len = 0;
	int g_len = 0;

	SHA256_CTX ctx;
	uint8_t buf[64] = {0};

	fn_t beta;
	fn_t gamma;
	fn_t alpha;
	fn_t xi;
	fn_t v;

	fn_t xi_w;
	fn_t v_pow;
	fn_t e1;
	fn_t e2;
	fn_t e3;

	int i;

	for (i = 0; i < 9; i++) {
		fn_set_zero(blind[i]);
	}
	fn_set_zero(fz[0]);
	fn_set_zero(fz[1]);
	fn_set_zero(beta);
	fn_set_zero(gamma);
	fn_set_zero(alpha);
	fn_set_zero(xi);
	fn_set_zero(v);
	fn_set_zero(xi_w);
	fn_set_zero(v_pow);
	fn_set_zero(e1);
	fn_set_zero(e2);
	fn_set_zero(e3);


	// Round 1
	/////////////////////////////////////////////////////////////////////////////
	printf("Round 1\n");

	// random blinding scalars
	printf("blinding scalars:\n");
	// constant_sage
	fn_from_hex(blind[0], sage_hex_blind0);
	fn_from_hex(blind[1], sage_hex_blind1);
	fn_from_hex(blind[2], sage_hex_blind2);
	fn_from_hex(blind[3], sage_hex_blind3);
	fn_from_hex(blind[4], sage_hex_blind4);
	fn_from_hex(blind[5], sage_hex_blind5);
	fn_from_hex(blind[6], sage_hex_blind6);
	fn_from_hex(blind[7], sage_hex_blind7);
	fn_from_hex(blind[8], sage_hex_blind8);
	for (i = 0; i < 9; i++) {
		//fn_rand(blind[i]);
		printf("  blind[%d] = ", i); fn_print("", blind[i]);
	}
	printf("\n");

	/*
	a(X) = interpolate(in[0..n-1])
	b(X) = interpolate(in[n..2n-1])
	c(X) = interpolate(in[2n..3n-1])

	NOTE: blindness order is different with the paper!
	a(X) += (blind[0] + blind[1] * X) * (X^n - 1)
	b(X) += (blind[2] + blind[3] * X) * (X^n - 1)
	c(X) += (blind[4] + blind[5] * X) * (X^n - 1)

	A = a(s) * G1
	B = b(s) * G1
	C = c(s) * G1
	*/
	for (i = 0; i < n; i++) {
		fn_copy(a[i], in[i]);
		fn_copy(b[i], in[n + i]);
		fn_copy(c[i], in[2 * n + i]);
	}
	poly_interpolate(a, H, nbits);
	poly_interpolate(b, H, nbits);
	poly_interpolate(c, H, nbits);
	poly_print(" a'(x) =", a, n);
	poly_print(" b'(x) =", b, n);
	poly_print(" c'(x) =", c, n);

	a_len = b_len = c_len = n;
	poly_add_blind(a, &a_len, blind, 2);
	poly_add_blind(b, &b_len, blind + 2, 2);
	poly_add_blind(c, &c_len, blind + 4, 2);
	poly_print(" a(x) =", a, a_len);
	poly_print(" b(x) =", b, b_len);
	poly_print(" c(x) =", c, c_len);

	point_multi_exponent(A, a, a_len, CRS);
	point_multi_exponent(B, b, a_len, CRS);
	point_multi_exponent(C, c, a_len, CRS);

	point_print(" A = a(s)G1 = ", A);
	point_print(" B = b(s)G1 = ", B);
	point_print(" C = c(s)G1 = ", C);
	printf("\n");

	// Round 2
	/////////////////////////////////////////////////////////////////////////////

	printf("Round 2\n");

	/*
	beta = Hash(A, B, C, H)
	gamma = Hash(A, B, C, H, beta)

	论文中的 (w_i)_{i in [l]} 是 [1, w, w^2, ...] ?
	*/
	/*
	sha256_init(&ctx);
	point_to_bytes(A, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(B, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(C, buf); sha256_update(&ctx, buf, 64);
	sha256_finish(&ctx, buf);
	fn_from_bytes(beta, buf);

	sha256_init(&ctx);
	point_to_bytes(A, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(B, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(C, buf); sha256_update(&ctx, buf, 64);
	fn_to_bytes(beta, buf); sha256_update(&ctx, buf, 32);
	sha256_finish(&ctx, buf);
	fn_from_bytes(gamma, buf);*/
	
	// constant_sage
	fn_from_hex(beta, sage_hex_beta);
	fn_from_hex(gamma, sage_hex_gamma);

	fn_print("beta = ", beta);
	fn_print("gamma = ", gamma);


	/*
	z[0] = 1
	z[1] = z[i-1]
		  beta * H[i]        + (a[i] + gamma)
		* -----------------------------------
		  beta * H[sigma(i)] + (a[i] + gamma)

		  beta * H[i] * k1   + (b[i] + gamma)
		* ---------------------------------------
		  beta * H[sigma(n + i)] + (b[i] + gamma)

		  beta * H[i] * k2   + (c[i] + gamma)
		* -----------------------------------------
		  beta * H[sigma(2*n + i)] + (c[i] + gamma)
	*/
	fn_set_one(z[0]);
	for (i = 1; i < n; i++) {
		int prod_is_one = 1;
		fn_t prod_num;
		fn_t prod_den;
		fn_t sum;
		fn_t num;
		fn_t den;

		if (sigma[i - 1] != i - 1) {
			fn_add(sum, in[i - 1], gamma);
			fn_mul(num, beta, H[i - 1]);
			fn_mul(den, beta, H[sigma[i - 1]]);
			fn_add(prod_num, num, sum);
			fn_add(prod_den, den, sum);
			prod_is_one = 0;
		}
		if (sigma[n + i - 1] != n + i - 1) {
			fn_add(sum, in[n + i - 1], gamma);
			fn_mul(num, beta, H[i - 1]);
			fn_mul(num, num, k1);
			fn_mul(den, beta, H[sigma[n + i - 1]]);
			fn_add(num, num, sum);
			fn_add(den, den, sum);
			if (prod_is_one) {
				fn_copy(prod_num, num);
				fn_copy(prod_den, den);
			} else {
				fn_mul(prod_num, prod_num, num);
				fn_mul(prod_den, prod_den, den);
			}
			prod_is_one = 0;
		}
		if (sigma[2*n + i - 1] != 2*n + i - 1) {
			fn_add(sum, in[2*n + i - 1], gamma);
			fn_mul(num, beta, H[i - 1]);
			fn_mul(num, num, k2);
			fn_mul(den, beta, H[sigma[2*n + i - 1]]);
			fn_add(num, num, sum);
			fn_add(den, den, sum);
			if (prod_is_one) {
				fn_copy(prod_num, num);
				fn_copy(prod_den, den);
			} else {
				fn_mul(prod_num, prod_num, num);
				fn_mul(prod_den, prod_den, den);
			}
			prod_is_one = 0;
		}

		if (prod_is_one) {
			fn_copy(z[i], z[i - 1]);
		} else {
			fn_inv(prod_den, prod_den);
			fn_mul(prod_num, prod_num, prod_den);
			fn_mul(z[i], z[i - 1], prod_num);
		}
	}

	/*
	z(X) = interpolate(z[0..n-1])
	z(X) += b[6] + b[7] * X + b[8] * X^2  # NOTE: blindness order is different with the paper!
	Z = z(s) * G1
	*/
	poly_interpolate(z, H, nbits);
	poly_print(" z'(x) =", z, n);

	z_len = n;
	poly_add_blind(z, &z_len, blind + 6, 3);
	poly_print(" z(x) =", z, z_len);

	point_multi_exponent(Z, z, z_len, CRS);
	point_print(" Z = z(s)*G1 = ", Z);


	// Round 3
	/////////////////////////////////////////////////////////////////////////////

	printf("\nRound 3\n");

	/*
	alpha = Hash(A, B, C, Z)
	*/
	/*
	sha256_init(&ctx);
	point_to_bytes(A, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(B, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(C, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(Z, buf); sha256_update(&ctx, buf, 64);
	sha256_finish(&ctx, buf);
	fn_from_bytes(alpha, buf);
	*/
	// constant_sage
	fn_from_hex(alpha, sage_hex_alpha);
	fn_print(" alpha = ", alpha);

	/*
	t(x) =   (qM(X) * b(x) + qL(X)) * a(x)
		+ qR(X) * b(x)
		+ qO(X) * c(x)
		+ PI(X)
		+ qC(X)
	order(t) = order(a) + order(b) + order(qM) = 2(n + 1) + (n - 1) = 3n + 1

	t(X) = t(X)/(X^n - 1)
	order(t) = 2n + 1
	*/
	poly_mul(t, &t_len, qM, n, b, b_len);
	poly_add(t, &t_len, t, t_len, qL, n);
	poly_mul(t, &t_len, t, t_len, a, a_len);
	poly_mul(r, &r_len, qR, n, b, b_len);
	poly_add(t, &t_len, t, t_len, r, r_len);
	poly_mul(r, &r_len, qO, n, c, c_len);
	poly_add(t, &t_len, t, t_len, r, r_len);
	//poly_add(t, &t_len, t, t_len, PI, n); // FIXME: how to handle PI ?
	poly_add(t, &t_len, t, t_len, qC, n);

	poly_div_ZH(t, &t_len, t, t_len, n);
	poly_print(" t(x) part-1", t, t_len);

	/*
	r(X)    = (a(X) + (gamma + beta * X))
		* (b(X) + (gamma + beta * k1 * X))
		* (c(X) + (gamma + beta * k2 * X))
		* z(X)
	order(r) = order(a) + order(b) + order(c) + order(z) = 3(n + 1) + (n + 2) = 4n + 5
	
	g(X)    = (S_sigma1(X) * beta + a(X) + gamma)
		* (S_sigma2(X) * beta + b(X) + gamma)
		* (S_sigma3(X) * beta + c(X) + gamma)
		* z(X * w)
	order(g) = 4n + 5

	z(X)   = a0 + a1*X   + a2*X^2 + ...
	z(X*w) = a0 + a1*X*w + a2*X^2*w^2 + ....

	f3(X) = ((r(X) - g(X)) * alpha)/(X^n - 1)
	order(f3) = 3n + 5

	t(X) = t(X) + f3(X)
	*/

	fn_copy(fz[0], gamma);
	fn_copy(fz[1], beta);
	poly_add(r, &r_len, a, a_len, fz, 2);
	fn_mul(fz[1], beta, k1);
	poly_add(f1, &f1_len, b, b_len, fz, 2);
	poly_mul(r, &r_len, r, r_len, f1, f1_len);
	fn_mul(fz[1], beta, k2);
	poly_add(f1, &f1_len, c, c_len, fz, 2);
	poly_mul(r, &r_len, r, r_len, f1, f1_len);
	poly_mul(r, &r_len, r, r_len, z, z_len);

	poly_mul_scalar(r, &r_len, r, r_len, alpha);
	
	// g
	
	poly_mul_scalar(g, &g_len, S_sigma1, n, beta);
	poly_add(g, &g_len, g, g_len, a, a_len);
	poly_add_scalar(g, &g_len, g, g_len, gamma);

	poly_mul_scalar(f1, &f1_len, S_sigma2, n, beta);
	poly_add(f1, &f1_len, f1, f1_len, b, b_len);
	poly_add_scalar(f1, &f1_len, f1, f1_len, gamma);
	poly_mul(g, &g_len, g, g_len, f1, f1_len);

	poly_mul_scalar(f1, &f1_len, S_sigma3, n, beta);
	poly_add(f1, &f1_len, f1, f1_len, c, c_len);
	poly_add_scalar(f1, &f1_len, f1, f1_len, gamma);
	poly_mul(g, &g_len, g, g_len, f1, f1_len);

	poly_mul_vector(f1, &f1_len, z, z_len, H, n);
	poly_mul(g, &g_len, g, g_len, f1, f1_len);

	poly_mul_scalar(g, &g_len, g, g_len, alpha);
	
	poly_sub(r, &r_len, r, r_len, g, g_len);
	
	
	poly_div_ZH(r, &r_len, r, r_len, n);
	poly_print(" t(x) part-3", r, r_len);

	poly_add(t, &t_len, t, t_len, r, r_len); // t = t + r

	/*
	r(X) = ((z(X) - 1) * L1(X) * alpha^2)/(X^n - 1)
	order(r) = (n + 2) + (n - 1) - n = n + 1

	t(X) = t(X) + r(X)
	*/
	fn_set_one(e1);
	poly_sub_scalar(r, &r_len, z, z_len, e1);
	poly_mul(r, &r_len, r, r_len, L1, n);
	fn_sqr(e1, alpha);
	poly_mul_scalar(r, &r_len, r, r_len, e1); // FIXME: 这里怎么没用L1(X)?
	poly_div_ZH(r, &r_len, r, r_len, n);
	poly_print(" t(x) part-4", r, r_len);

	poly_add(t, &t_len, t, t_len, r, r_len);
	poly_print(" t(x)", t, t_len);

	/*
	t(X) = t_lo(X) + t_mid(X) * X^n + t_hi(X) * X^2n
	order(t) = 3n + 5

	T_lo = t_lo(s) * G1
	T_mid = t_mid(s) * G1
	T_hi = t_hi(s) * G1
	*/
	point_multi_exponent(T_lo, t, n, CRS);
	point_multi_exponent(T_mid, t + n, n, CRS);
	point_multi_exponent(T_hi, t + n + n, t_len - n - n, CRS);
	point_print(" T_lo  = t_lo(s)G1  = ", T_lo);
	point_print(" T_mid = t_mid(s)G1 = ", T_mid);
	point_print(" T_hi  = t_hi(s)G1  = ", T_hi);


	// Round 4
	/////////////////////////////////////////////////////////////////////////////

	printf("\nRound 4\n");

	// xi = Hash(A, B, C, Z, T_lo, T_mid, T_hi)
	/*
	sha256_init(&ctx);
	point_to_bytes(A, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(B, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(C, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(Z, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(T_lo, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(T_mid, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(T_hi, buf); sha256_update(&ctx, buf, 64);
	sha256_finish(&ctx, buf);
	fn_from_bytes(xi, buf);
	*/
	// constant_sage
	fn_from_hex(xi, sage_hex_xi);
	fn_print(" xi        = ", xi);


	/*
	_a = a(xi)
	_b = b(xi)
	_c = c(xi)
	_s_sigma1 = S_sigma1(xi)
	_s_sigma2 = S_sigma2(xi)
	_t = t(xi)
	_z_w = z(xi * w) = z_w(xi)
	*/
	poly_eval(_a, a, a_len, xi);
	poly_eval(_b, b, b_len, xi);
	poly_eval(_c, c, c_len, xi);
	poly_eval(_s_sigma1, S_sigma1, n, xi);
	poly_eval(_s_sigma2, S_sigma2, n, xi);
	poly_eval(_t, t, t_len, xi);
	fn_mul(xi_w, xi, H[1]); // H[1] == w
	poly_eval(_z_w, z, z_len, xi_w);
	fn_print(" _a        = ", _a);
	fn_print(" _b        = ", _b);
	fn_print(" _c        = ", _c);
	fn_print(" _s_sigma1 = ", _s_sigma1);
	fn_print(" _s_sigma2 = ", _s_sigma2);
	fn_print(" _t        = ", _t);
	fn_print(" _z_w      = ", _z_w);
	printf("\n");

	printf(" r(x)\n");

	/*
	e1 = (beta * xi      + _a + gamma)
		* (beta * xi * k1 + _b + gamma)
		* (beta * xi * k2 + _c + gamma)
	*/
	fn_mul(e2, beta, xi);
	fn_add(e1, e2, _a);
	fn_add(e1, e1, gamma);
	fn_mul(e3, e2, k1);
	fn_add(e3, e3, _b);
	fn_add(e3, e3, gamma);
	fn_mul(e1, e1, e3);
	fn_mul(e2, e2, k2);
	fn_add(e2, e2, _c);
	fn_add(e2, e2, gamma);
	fn_mul(e1, e1, e2);
	fn_print("  e1 = ", e1);

	/*
	e2 =      (beta * _s_sigma1 + _a + gamma)
		* (beta * _s_sigma2 + _b + gamma)
		* beta * _z_w

	e3 = L_1(xi) * alpha
	*/
	fn_mul(e2, beta, _s_sigma1);
	fn_add(e2, e2, _a);
	fn_add(e2, e2, gamma);
	fn_mul(e3, beta, _s_sigma2);
	fn_add(e3, e3, _b);
	fn_add(e3, e3, gamma);
	fn_mul(e2, e2, e3);
	fn_mul(e2, e2, beta);
	fn_mul(e2, e2, _z_w);
	fn_print("  e2 = ", e2);

	poly_eval_L1(e3, n, xi);
	fn_mul(e3, e3, alpha);
	fn_print("  e3 = ", e3);

	/*
	r(X) = (z(X) * (e1 + e3) - S_sigma3(X) * e2) * alpha
		----------------   ----------------
		   n + 3                   n

	f1(X) = qM(X) * (_a * _b)
		+ qL(X) * _a
		+ qR(X) * _b
		+ qO(X) * _c
		+ qC(X)

	r(X) = f1(X) + r(X)
	*/
	fn_add(e1, e1, e3);
	poly_mul_scalar(r, &r_len, z, z_len, e1);
	poly_mul_scalar(f1, &f1_len, S_sigma3, n, e2);
	poly_sub(r, &r_len, r, r_len, f1, f1_len);
	poly_mul_scalar(r, &r_len, r, r_len, alpha);
	poly_print(" r(x) = ", r, r_len);

	fn_mul(e1, _a, _b);
	poly_mul_scalar(f1, &f1_len, qM, n, e1);
	poly_mul_scalar(f2, &f2_len, qL, n, _a);
	poly_add(f1, &f1_len, f1, f1_len, f2, f2_len);
	poly_mul_scalar(f2, &f2_len, qR, n, _b);
	poly_add(f1, &f1_len, f1, f1_len, f2, f2_len);
	poly_mul_scalar(f2, &f2_len, qO, n, _c);
	poly_add(f1, &f1_len, f1, f1_len, f2, f2_len);
	poly_add(f1, &f1_len, f1, f1_len, qC, n);
	poly_print(" f1(x) = ", f1, f1_len);

	poly_add(r, &r_len, f1, f1_len, r, r_len);
	poly_print(" final r(x) = ", r, r_len);

	poly_eval(_r, r, r_len, xi);
	fn_print(" _r = ", _r);


	// Round 5:
	/////////////////////////////////////////////////////////////////////////////

	printf("\nRound 5\n");

	/*
	v = Hash(A, B, C, Z, T,
		_a, _b, _c, _s_sigma1, _s_sigma2, _z_w,
		T_lo, T_mid, T_hi, _r)
	*/
	/*
	sha256_init(&ctx);
	point_to_bytes(A, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(B, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(C, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(Z, buf); sha256_update(&ctx, buf, 64);
	fn_to_bytes(_a, buf); sha256_update(&ctx, buf, 32);
	fn_to_bytes(_b, buf); sha256_update(&ctx, buf, 32);
	fn_to_bytes(_c, buf); sha256_update(&ctx, buf, 32);
	fn_to_bytes(_s_sigma1, buf); sha256_update(&ctx, buf, 32);
	fn_to_bytes(_s_sigma2, buf); sha256_update(&ctx, buf, 32);
	fn_to_bytes(_z_w, buf); sha256_update(&ctx, buf, 32);
	point_to_bytes(T_lo, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(T_mid, buf); sha256_update(&ctx, buf, 64);
	point_to_bytes(T_hi, buf); sha256_update(&ctx, buf, 64);
	fn_to_bytes(_r, buf); sha256_update(&ctx, buf, 32);
	sha256_finish(&ctx, buf);
	fn_from_bytes(v, buf);
	*/
	// constant_sage
	fn_from_hex(v, sage_hex_v);
	fn_print(" v = ", v);


	/*
	t2(X) = t0(X) + t1(X) * xi^n + t2(X) * (xi^n)^2 - _t	# = n + 5
		+ (r(X) - _r) * v				# = n + 2
		+ (a(X) - _a) * v^2				# = n + 1
		+ (b(X) - _b) * v^3				# = n + 1
		+ (c(X) - _c) * v^4				# = n + 1
		+ (S_sigma1(X) - _s_sigma1) * v^5		# = n - 1
		+ (S_sigma2(X) - _s_sigma2) * v^6		# = n - 1

	t2(X) = t2(X)/(X - xi) = W_xi(X)			# = n + 4

	W_xi = t2(X) * CRS = [W_xi]_1
	*/

	fn_t *t0 = t;
	fn_t *t1 = t + n;
	fn_t *t2 = t + n + n;
	int t0_len = n;
	int t1_len = n;
	int t2_len = t_len - n - n;


	fn_exp(e1, xi, (uint32_t *)&n, 1);
	poly_mul_scalar(t1, &t1_len, t1, t1_len, e1);
	fn_sqr(e1, e1);
	poly_mul_scalar(t2, &t2_len, t2, t2_len, e1);
	poly_add(t1, &t1_len, t1, t1_len, t0, t0_len);
	poly_add(t2, &t2_len, t2, t2_len, t1, t1_len);
	poly_sub_scalar(t2, &t2_len, t2, t2_len, _t);
	poly_print(" W_xi(X) step-1", t2, t2_len);

	fn_copy(v_pow, v);
	poly_sub_scalar(r, &r_len, r, r_len, _r);
	poly_mul_scalar(r, &r_len, r, r_len, v_pow);
	poly_add(t2, &t2_len, t2, t2_len, r, r_len);
	poly_print(" W_xi(X) step-2", t2, t2_len);

	fn_sqr(v_pow, v);
	poly_sub_scalar(a, &a_len, a, a_len, _a);
	poly_mul_scalar(a, &a_len, a, a_len, v_pow);
	poly_add(t2, &t2_len, t2, t2_len, a, a_len);
	poly_print(" W_xi(X) step-3", t2, t2_len);

	fn_mul(v_pow, v_pow, v);
	poly_sub_scalar(b, &b_len, b, b_len, _b);
	poly_mul_scalar(b, &b_len, b, b_len, v_pow);
	poly_add(t2, &t2_len, t2, t2_len, b, b_len);
	poly_print(" W_xi(X) step-4", t2, t2_len);

	fn_mul(v_pow, v_pow, v);
	poly_sub_scalar(c, &c_len, c, c_len, _c);
	poly_mul_scalar(c, &c_len, c, c_len, v_pow);
	poly_add(t2, &t2_len, t2, t2_len, c, c_len);
	poly_print(" W_xi(X) step-5", t2, t2_len);

	fn_mul(v_pow, v_pow, v);
	poly_sub_scalar(t0, &t0_len, S_sigma1, n, _s_sigma1);
	poly_mul_scalar(t0, &t0_len, t0, t0_len, v_pow);
	poly_add(t2, &t2_len, t2, t2_len, t0, t0_len);
	poly_print(" W_xi(X) step-6", t2, t2_len);

	fn_mul(v_pow, v_pow, v);
	poly_sub_scalar(t1, &t1_len, S_sigma2, n, _s_sigma2);
	poly_mul_scalar(t1, &t1_len, t1, t1_len, v_pow);
	poly_add(t2, &t2_len, t2, t2_len, t1, t1_len);
	poly_print(" W_xi(X) step-7", t2, t2_len);

	poly_div_x_sub_scalar(t2, &t2_len, xi);
	poly_print(" W_xi(X) step-8", t2, t2_len);

	point_multi_exponent(W_xi, t2, t2_len, CRS);
	point_print(" W_xi = W_xi(s)*G1     = ", W_xi);

	/*
	z(x) : n + 3
	                 (z(X) - _z_w)
	w_{xi w}(X) = -------------
	                 X - xi * w

	W_detla_w = w_{xi w} * H
	*/
	poly_sub_scalar(z, &z_len, z, z_len, _z_w);
	poly_div_x_sub_scalar(z, &z_len, xi_w);
	point_multi_exponent(W_xi_w, z, z_len, CRS); // z_len == 2
	point_print(" W_xi_w = W_xi_w(s)*G1 = ", W_xi_w);

	printf("\n");
	printf("proof:\n");
	point_print(" A         = ", A);
	point_print(" B         = ", B);
	point_print(" C         = ", C);
	point_print(" Z         = ", Z);
	point_print(" T_lo      = ", T_lo);
	point_print(" T_mid     = ", T_mid);
	point_print(" T_hi      = ", T_hi);
	point_print(" W_xi      = ", W_xi);
	point_print(" W_xi_w    = ", W_xi_w);
	fn_print(   " _a        = ", _a);
	fn_print(   " _b        = ", _b);
	fn_print(   " _c        = ", _c);
	fn_print(   " _s_sigma1 = ", _s_sigma1);
	fn_print(   " _s_sigma2 = ", _s_sigma2);
	fn_print(   " _r        = ", _r);
	fn_print(   " _z_w      = ", _z_w);

	// FIXME: 输入错误数据，无法证明，返回 -1
	return 1;
}

