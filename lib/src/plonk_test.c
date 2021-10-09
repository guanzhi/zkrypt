#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bn254.h"
#include "bn254_test.h"
#include "plonk.h"


int test_plonk(void)
{
	int n = 16;
	fn_t s;
	affine_point_t *CRS = NULL;
	fn_t *H = poly_new(16);
	fn_t k1;
	fn_t k2;

	fn_t *L1 = poly_new(16);
	fn_t *qM = poly_new(16);
	fn_t *qL = poly_new(16);
	fn_t *qR = poly_new(16);
	fn_t *qO = poly_new(16);
	fn_t *qC = poly_new(16);
	fn_t *PI = poly_new(16);
	int sigma[16 * 3];
	fn_t *S_sigma1 = poly_new(16);
	fn_t *S_sigma2 = poly_new(16);
	fn_t *S_sigma3 = poly_new(16);
	fn_t *in = poly_new(16 * 3);
	point_t A;
	point_t B;
	point_t C;
	point_t Z;
	point_t T_lo;
	point_t T_mid;
	point_t T_hi;
	point_t W_delta;
	point_t W_delta_w;
	fn_t _a;
	fn_t _b;
	fn_t _c;
	fn_t _s_sigma1;
	fn_t _s_sigma2;
	fn_t _t;
	fn_t _r;
	fn_t _z_w;

	int i;
	printf("%s\n", __FUNCTION__);

	fn_rand(s);
	CRS = (affine_point_t *)malloc(sizeof(affine_point_t) * (n + 6));
	plonk_generate_crs(CRS, s, n + 6);

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

	fn_set_word(k1, 2);
	fn_set_word(k2, 3);

	poly_rand(L1, 16);
	poly_rand(qM, 16);
	poly_rand(qL, 16);
	poly_rand(qR, 16);
	poly_rand(qO, 16);
	poly_rand(qC, 16);
	poly_rand(PI, 16);

	for (i = 0; i < 16; i++) {
		sigma[i] = 15 - i;
	}
	poly_rand(S_sigma1, 16);
	poly_rand(S_sigma2, 16);
	poly_rand(in, 16 * 3);


	plonk_proof(CRS, H, 4, k1, k2, L1, qM, qL, qR, qO, qC, PI,
		sigma, S_sigma1, S_sigma2, S_sigma3, in,
		&A, &B, &C, &Z, &T_lo, &T_mid, &T_hi, &W_delta, &W_delta_w,
		_a, _b, _c, _s_sigma1, _s_sigma2, _t, _r, _z_w);





	return 0;
}

int test_plonk_from_sage_simple(void)
{
	int n = 4;
	fn_t s;
	affine_point_t *CRS = NULL;
	fn_t *H = poly_new(n * 3);
	fn_t k1;
	fn_t k2;

	fn_t *L1 = poly_new(n);
	fn_t *qM = poly_new(n);
	fn_t *qL = poly_new(n);
	fn_t *qR = poly_new(n);
	fn_t *qO = poly_new(n);
	fn_t *qC = poly_new(n);
	fn_t *PI = poly_new(n);
	int sigma[4 * 3] = {4, 5, 6, 8, 0, 1, 2, 9, 3, 7, 11, 10};
	fn_t *S_sigma1 = poly_new(n);
	fn_t *S_sigma2 = poly_new(n);
	fn_t *S_sigma3 = poly_new(n);
	fn_t *in = poly_new(n * 3);
	point_t A;
	point_t B;
	point_t C;
	point_t Z;
	point_t T_lo;
	point_t T_mid;
	point_t T_hi;
	point_t W_delta;
	point_t W_delta_w;
	fn_t _a;
	fn_t _b;
	fn_t _c;
	fn_t _s_sigma1;
	fn_t _s_sigma2;
	fn_t _t;
	fn_t _r;
	fn_t _z_w;

	int i;
	printf("%s\n", __FUNCTION__);

	fn_from_hex(s, sage_hex_s);
	
	CRS = (affine_point_t *)malloc(sizeof(affine_point_t) * (n + 6));
	plonk_generate_crs(CRS, s, n + 6);

	fn_set_word(H[0], 1);
	fn_from_hex(H[1], sage_hex_w0);
	for(i=2;i<n;++i) {
		fn_mul(H[i], H[i-1], H[1]);
		//fn_print("H ", H[i]);
	}

	fn_set_word(k1, 5);
	fn_set_word(k2, 7);
	
	for(i=0;i<n;++i) {
		fn_mul(H[n + i], H[i], k1);
		fn_mul(H[2*n + i], H[i], k2);
	}
	
	int vqL[4] = {0, 0, 0, 1};
	
	for(i=0;i<4;++i) {
		fn_set_word(qL[i], vqL[i]);
		fn_set_word(qR[i], vqL[i]);
		fn_from_hex(qO[i], sage_hex_neg_1);
		fn_set_word(qM[i], 1 - vqL[i]);
		fn_set_word(qC[i], 0);
		fn_set_word(PI[i], 0);
	}
	poly_interpolate(qL, H, 2);
	poly_interpolate(qR, H, 2);
	poly_interpolate(qO, H, 2);
	poly_interpolate(qM, H, 2);
	poly_interpolate(qC, H, 2);
	poly_interpolate(PI, H, 2);
	
	fn_set_word(L1[0], 1);
	fn_set_word(L1[1], 0);
	fn_set_word(L1[2], 0);
	fn_set_word(L1[3], 0);
	poly_interpolate(L1, H, 2);
	
	for(i=0;i<4;++i) {
		fn_copy(S_sigma1[i], H[sigma[i]]);
		fn_copy(S_sigma2[i], H[sigma[n+i]]);
		fn_copy(S_sigma3[i], H[sigma[2*n+i]]);
	}
	poly_interpolate(S_sigma1, H, 2);
	poly_interpolate(S_sigma2, H, 2);
	poly_interpolate(S_sigma3, H, 2);
	
	int vin[12] = {3, 4, 5, 9, 3, 4, 5, 16, 9, 16, 25, 25};
	// int vin[12] = {12, 5, 13, 144, 12, 5, 13, 25, 144, 25, 169, 169};
	for(i=0;i<12;++i) {
		fn_set_word(in[i], vin[i]);
	}

	plonk_proof(CRS, H, 2, k1, k2, L1, qM, qL, qR, qO, qC, PI,
		sigma, S_sigma1, S_sigma2, S_sigma3, in,
		&A, &B, &C, &Z, &T_lo, &T_mid, &T_hi, &W_delta, &W_delta_w,
		_a, _b, _c, _s_sigma1, _s_sigma2, _t, _r, _z_w);
	
	return 0;
	/* Answer:
	015076413a1e3ccd8f8e06a3e3f3593473c7fb6c75ebee4f481465cb91f0a858 189ff80237a68c56e9a04255b09b36bfff44d0643b2551aa8a7fe7f50319bc0a
	258f280f6c3b4d89047293b190884e674350c8b5105a95b35e8500d57fb66dee 1100e2e49a1a696d0dbad87b07763f631abb256ec7613ef3112066f193b6bf93
	209cb88857f7f3691b389559bbd92e46c73f69a96ce30fed4ded56e988f2777c 10786bdf254bb11422e65a019ef465f660f00e7bb815c45c351b3bdc0b0430dc
	056e290547f77e93184a191aea1dd4f7e8b9a6628caa3478bf93a0c138cfb72d 11b027fba02dc51a7f8d092b3d2e9c32f205538abbce60ffadbfd9a7ef6d5577
	0b077787871b2a77242c4064854fcbedb2f0397ef9525f0d92afac3539dce2d5 0d99aeb6fd0f712891ec7a01f27dd514ca4a194bf93655ee6f53064747e5c83f
	16ff27db58ee317a728125be457db2ae65da591663211bfed0146b7114b11c61 2bd82bfa91d559705b1fc4153bd149524106fb244fde83f444a50b730c1f1b5e
	0e624a95aaa9973d0d04a41abfcd187944c65a10ef3f43f0574b05fe56c30414 08a1c1db8b10604f7980d64ad121d8b806bd407183bae466eb5769a3ccbbf519
	0631d246d568f2dcda81257288be0abe5b6be4ec33b16383bdfd945f3dbf8147 0c812cd688f72b74ca5b79fcdf5bcf146f2da816134bc77be5e64de69e21b1ab
	09d8cad01d18aae6a141abea7f4bf0d05a9f1dfca12d38f7766c53189b5f42b4 0ae73daef13721014524ac695faefdc4d523849f1391af236a87505fa0f66e17
	21b575ce967712d9add2e169321307b03f6a99fb0beebd7dd1a63402833f5981
	02847226776e91d204387f0d7068644b77c1e12c751907c6ffe32e54067cf549
	013dae72fb447023df84e4bd41fbfe7efaa26fb9fc3feedd53186f58c9973c92
	125e521ea6c2acab9d00b1ae4bdd177ff39c5b2352e572af392ad3869ba39ed0
	1efb4c64bf9f43d86b583cd3ad3020718fe89c9b11aa726456dae6671eaea387
	0f68d978fb6fb49b3dac4e86d2638e571cb1a5ced0e8ab9b9d69dd25a5868ff7
	06421e1705bbbeaa3098d0afe2b85b3bc0443d1de9493b7a7b1b53b0d375490a
	*/
}

int test_plonk_from_sage_mimc_small(void)
{
	int round = 12; // setup
	int n = 64; // setup
	fn_t s;
	affine_point_t *CRS = NULL;
	fn_t *H = poly_new(n * 3);
	fn_t k1;
	fn_t k2;

	fn_t *L1 = poly_new(n);
	fn_t *qM = poly_new(n);
	fn_t *qL = poly_new(n);
	fn_t *qR = poly_new(n);
	fn_t *qO = poly_new(n);
	fn_t *qC = poly_new(n);
	fn_t *PI = poly_new(n);
	fn_t *wit = poly_new(5*round+3);
	fn_t *constant = poly_new(round);
	int sigma[64*3] = {}; // setup
	int pos[63] = {}, start[63] = {};  // setup: 5 * round + 3
	fn_t *S_sigma1 = poly_new(n);
	fn_t *S_sigma2 = poly_new(n);
	fn_t *S_sigma3 = poly_new(n);
	fn_t *in = poly_new(n * 3);
	point_t A;
	point_t B;
	point_t C;
	point_t Z;
	point_t T_lo;
	point_t T_mid;
	point_t T_hi;
	point_t W_delta;
	point_t W_delta_w;
	fn_t _a;
	fn_t _b;
	fn_t _c;
	fn_t _s_sigma1;
	fn_t _s_sigma2;
	fn_t _t;
	fn_t _r;
	fn_t _z_w;

	int i, j;
	printf("%s\n", __FUNCTION__);

	fn_from_hex(s, sage_hex_s);
	
	CRS = (affine_point_t *)malloc(sizeof(affine_point_t) * (n + 6));
	plonk_generate_crs(CRS, s, n + 6);

	fn_set_word(H[0], 1);
	fn_from_hex(H[1], sage_hex_w_mimc_small); // setup
	for(i=2;i<n;++i) {
		fn_mul(H[i], H[i-1], H[1]);
		//fn_print("H ", H[i]);
	}

	fn_set_word(k1, 5); // setup
	fn_set_word(k2, 7); // setup
	
	for(i=0;i<n;++i) {
		fn_mul(H[n + i], H[i], k1);
		fn_mul(H[2*n + i], H[i], k2);
	}
	
	fn_from_hex(constant[0], sage_hex_c0); // setup
	fn_from_hex(constant[1], sage_hex_c1);
	fn_from_hex(constant[2], sage_hex_c2);
	fn_from_hex(constant[3], sage_hex_c3);
	fn_from_hex(constant[4], sage_hex_c4);
	fn_from_hex(constant[5], sage_hex_c5);
	fn_from_hex(constant[6], sage_hex_c6);
	fn_from_hex(constant[7], sage_hex_c7);
	fn_from_hex(constant[8], sage_hex_c8);
	fn_from_hex(constant[9], sage_hex_c9);
	fn_from_hex(constant[10], sage_hex_c10);
	fn_from_hex(constant[11], sage_hex_c11);
	
	int vqL[4] = {1, 0, 0, 1}; // addition gate = 1, multi gate = 0
	int va[64] = {}, vb[64] = {}, vc[64] = {}; // setup
	
	for(j=0;j<round;++j) {
		for(i=0;i<4;++i) {
			fn_set_word(qL[5*j+i], vqL[i]);
			fn_set_word(qR[5*j+i], vqL[i]);
			fn_from_hex(qO[5*j+i], sage_hex_neg_1);
			fn_set_word(qM[5*j+i], 1 - vqL[i]);
			fn_set_word(qC[5*j+i], 0);
			fn_set_word(PI[5*j+i], 0);
		}
		fn_from_hex(qL[5*j+4], sage_hex_neg_1);
		fn_set_word(qR[5*j+4], 0);
		fn_set_word(qO[5*j+4], 0);
		fn_set_word(qM[5*j+4], 0);
		fn_copy(qC[5*j+4], constant[j]);
		fn_set_word(PI[5*j+4], 0);
	}
	
	for(i=5*round;i<n;++i) {
		fn_set_word(qL[i], 0);
		fn_set_word(qR[i], 0);
		fn_set_word(qO[i], 0);
		fn_set_word(qM[i], 0);
		fn_set_word(qC[i], 0);
		fn_set_word(PI[i], 0);
		va[i] = 0;
		vb[i] = 0;
		vc[i] = 0;
	}
	
	for(j=0;j<round;++j) {
		va[5*j] = 5*j+2;
		va[5*j+1] = 5*j+4;
		va[5*j+2] = 5*j+4;
		va[5*j+3] = 5*j+6;
		va[5*j+4] = 5*j+3;
		vb[5*j] = 5*j+3;
		vb[5*j+1] = 5*j+4;
		vb[5*j+2] = 5*j+5;
		vb[5*j+3] = 5*j-3;
		vb[5*j+4] = 0;
		vc[5*j] = 5*j+4;
		vc[5*j+1] = 5*j+5;
		vc[5*j+2] = 5*j+6;
		vc[5*j+3] = 5*j+7;
		vc[5*j+4] = 0;
	}
	vb[3] = 1;
		
	poly_interpolate(qL, H, 6); // setup
	poly_interpolate(qR, H, 6);
	poly_interpolate(qO, H, 6);
	poly_interpolate(qM, H, 6);
	poly_interpolate(qC, H, 6);
	poly_interpolate(PI, H, 6);
	
	
	memset(sigma, -1, sizeof(sigma));
	memset(pos, -1, sizeof(pos));
	memset(start, -1, sizeof(start));
	
	for(i=0;i<n;++i) {
		if(pos[va[i]] == -1) start[va[i]] = i;
		else sigma[i] = pos[va[i]];
		pos[va[i]] = i;
	}
	
	for(i=n;i<2*n;++i) {
		if(pos[vb[i-n]] == -1) start[vb[i-n]] = i;
		else sigma[i] = pos[vb[i-n]];
		pos[vb[i-n]] = i;
	}
	
	for(i=2*n;i<3*n;++i) {
		if(pos[vc[i-2*n]] == -1) start[vc[i-2*n]] = i;
		else sigma[i] = pos[vc[i-2*n]];
		pos[vc[i-2*n]] = i;
	}
	
	for(i=0;i<5*round+3;++i) {
		if(start[i] != -1) sigma[start[i]] = pos[i];
	}
	
	fn_set_word(L1[0], 1);
	for(i=1;i<n;++i) {
		fn_set_word(L1[i], 0);
	}
	poly_interpolate(L1, H, 6);
	
	for(i=0;i<n;++i) {
		fn_copy(S_sigma1[i], H[sigma[i]]);
		fn_copy(S_sigma2[i], H[sigma[n+i]]);
		fn_copy(S_sigma3[i], H[sigma[2*n+i]]);
	}
	poly_interpolate(S_sigma1, H, 6); // setup
	poly_interpolate(S_sigma2, H, 6);
	poly_interpolate(S_sigma3, H, 6);
	
	fn_t xl;
	fn_t xr;
	
	fn_from_hex(xr, sage_hex_xr);
	fn_from_hex(xl, sage_hex_xl);
	fn_set_word(wit[0], 0);
	fn_copy(wit[1], xr);
	fn_copy(wit[2], xl);
	
	for(j=0;j<round;++j) {
		fn_copy(wit[5*j+3], constant[j]);
		fn_add(wit[5*j+4], xl, constant[j]);
		fn_sqr(wit[5*j+5], wit[5*j+4]);
		fn_mul(wit[5*j+6], wit[5*j+5], wit[5*j+4]);
		fn_add(wit[5*j+7], wit[5*j+6], xr);
		fn_copy(xr, xl);
		fn_copy(xl, wit[5*j+7]);
	}
	
	for(i=0;i<n;++i) {
		fn_copy(in[i], wit[va[i]]);
		fn_copy(in[n+i], wit[vb[i]]);
		fn_copy(in[2*n+i], wit[vc[i]]);
	}

	plonk_proof(CRS, H, 6, k1, k2, L1, qM, qL, qR, qO, qC, PI,
		sigma, S_sigma1, S_sigma2, S_sigma3, in,
		&A, &B, &C, &Z, &T_lo, &T_mid, &T_hi, &W_delta, &W_delta_w,
		_a, _b, _c, _s_sigma1, _s_sigma2, _t, _r, _z_w);
	
	return 0;
	/* Answer:
	138aa3a294293c0f031043598f6e0852aeac796c15d1b8f4bda303f2d97bc1d8 018f308505e0bba2a893e04c8db23d79b8090f5e1629e4e052a3a1f0652b087a
	0177dd7c0c65d3375d4cd9079160deaa51bef4f0d4ce29d24bfe41eeef1a8fd2 2a2bbd0fcbeefecd736583f00159c6fefc0dbd2a5ec6fe1778abd43ae8d33ee4
	0f3f1d23410f8e0971cbd83e2efd163f5ce144dcfde75c7fbd3c244a5ced0073 17bb3c4de0e6f06dfd700b851fa58effed828c93a4217ba8c44984261ef58f24
	22b98324a8711ece6ceda494db1894d4e744b895ae802ec37bc6ff6bdc9ca660 2568b475c00096f138ba1a1f5cad703b911b6b53c9374ffcd1cfc619ac6fcb58
	1f26d73c7a43f7936555a2d9bc348fc74ed1029d8a7c9260ba071d3e6877dd66 0e5dee0095f3f8d629e9f20c2fc8e23a47b19c8392e787f673eca5b9f207db90
	20315a7624f7d459bd20d3f495ccab005e5e7a220ed8248397a4dc0446dc6e98 2186c64e716e7915a04d4aa9f53ba76f53153ed07afe0ddb85a0f7df66e92ef7
	20f2c2c319f8c69cf9b61114d47b9c0a2ce539cc646e5bf69a87cfc4c435fcab 2a5e64d7f97de26990d7627a3ec71902491e17641ee977c225c3671bccbe39c7
	05c75af734fd63767fc9e0e3174c1ab1bcc03ff4bb2af36e415facbac92219b3 1f57aa4638081adc2e6897912a31724440a618aae22b225355b659f57a356450
	21de8ad2fee39e14dcb9690b7d6210b34bafaf8fcb282e68bc71ca62f72ee673 06bb8194b49ec0cc540dc0c3fa32e34c19b964c106778d33882614642fe21506
	0c097d0f5fdaba6cbbc6f63f0fa9364de58256c6a13d679d01ec3cd567cd551a
	1ebbe54c842a7da33ed0fe39e0c2ceb67f20eada57b250d1f194e78cdffee2cd
	2c68111723bfd55da463de170485ccd906e9030258dd3f959d1e96521cde3618
	16a00d33ef81bd7f062153301a6633591077b4fe83d1807efc24eafdf31b0215
	1cb352adaecd3e567f5d2049f2f30e83a34b5184f0124991d5dd6602d413839e
	14417c61b2e0dcc35bbe43a91f2fa41e9b4d142e401eb1a676a0295a5f38fbf6
	1accc7b64339b0ce519a663543900d309ff579f15e4d171be7962051408f4053
	*/
}

int test_plonk_from_sage_mimc_medium(void)
{
	int round = 800; // setup
	int n = 4096; // setup
	fn_t s;
	affine_point_t *CRS = NULL;
	fn_t *H = poly_new(n * 3);
	fn_t k1;
	fn_t k2;

	fn_t *L1 = poly_new(n);
	fn_t *qM = poly_new(n);
	fn_t *qL = poly_new(n);
	fn_t *qR = poly_new(n);
	fn_t *qO = poly_new(n);
	fn_t *qC = poly_new(n);
	fn_t *PI = poly_new(n);
	fn_t *wit = poly_new(5 * round + 3);
	fn_t *constant = poly_new(round);
	int sigma[4096 * 3] = {}; // setup
	int pos[4003] = {}, start[4003] = {};  // setup: 5 * round + 3
	fn_t *S_sigma1 = poly_new(n);
	fn_t *S_sigma2 = poly_new(n);
	fn_t *S_sigma3 = poly_new(n);
	fn_t *in = poly_new(n * 3);
	point_t A;
	point_t B;
	point_t C;
	point_t Z;
	point_t T_lo;
	point_t T_mid;
	point_t T_hi;
	point_t W_delta;
	point_t W_delta_w;
	fn_t _a;
	fn_t _b;
	fn_t _c;
	fn_t _s_sigma1;
	fn_t _s_sigma2;
	fn_t _t;
	fn_t _r;
	fn_t _z_w;

	int i, j;
	printf("%s\n", __FUNCTION__);

	fn_from_hex(s, sage_hex_s);
	
	CRS = (affine_point_t *)malloc(sizeof(affine_point_t) * (n + 6));
	plonk_generate_crs(CRS, s, n + 6);

	fn_set_word(H[0], 1);
	fn_from_hex(H[1], sage_hex_w_mimc_medium); // setup
	for(i=2;i<n;++i) {
		fn_mul(H[i], H[i-1], H[1]);
		//fn_print("H ", H[i]);
	}

	fn_set_word(k1, 5); // setup
	fn_set_word(k2, 7); // setup
	
	for(i=0;i<n;++i) {
		fn_mul(H[n + i], H[i], k1);
		fn_mul(H[2*n + i], H[i], k2);
	}
	
	for(i=0;i<round;++i) fn_rand(constant[i]);
	
	int vqL[4] = {1, 0, 0, 1}; // addition gate = 1, multi gate = 0
	int va[4096] = {}, vb[4096] = {}, vc[4096] = {}; // setup
	
	for(j=0;j<round;++j) {
		for(i=0;i<4;++i) {
			fn_set_word(qL[5*j+i], vqL[i]);
			fn_set_word(qR[5*j+i], vqL[i]);
			fn_from_hex(qO[5*j+i], sage_hex_neg_1);
			fn_set_word(qM[5*j+i], 1 - vqL[i]);
			fn_set_word(qC[5*j+i], 0);
			fn_set_word(PI[5*j+i], 0);
		}
		fn_from_hex(qL[5*j+4], sage_hex_neg_1);
		fn_set_word(qR[5*j+4], 0);
		fn_set_word(qO[5*j+4], 0);
		fn_set_word(qM[5*j+4], 0);
		fn_copy(qC[5*j+4], constant[j]);
		fn_set_word(PI[5*j+4], 0);
	}
	
	for(i=5*round;i<n;++i) {
		fn_set_word(qL[i], 0);
		fn_set_word(qR[i], 0);
		fn_set_word(qO[i], 0);
		fn_set_word(qM[i], 0);
		fn_set_word(qC[i], 0);
		fn_set_word(PI[i], 0);
		va[i] = 0;
		vb[i] = 0;
		vc[i] = 0;
	}
	
	for(j=0;j<round;++j) {
		va[5*j] = 5*j+2;
		va[5*j+1] = 5*j+4;
		va[5*j+2] = 5*j+4;
		va[5*j+3] = 5*j+6;
		va[5*j+4] = 5*j+3;
		vb[5*j] = 5*j+3;
		vb[5*j+1] = 5*j+4;
		vb[5*j+2] = 5*j+5;
		vb[5*j+3] = 5*j-3;
		vb[5*j+4] = 0;
		vc[5*j] = 5*j+4;
		vc[5*j+1] = 5*j+5;
		vc[5*j+2] = 5*j+6;
		vc[5*j+3] = 5*j+7;
		vc[5*j+4] = 0;
	}
	vb[3] = 1;
		
	poly_interpolate(qL, H, 12); // setup
	poly_interpolate(qR, H, 12);
	poly_interpolate(qO, H, 12);
	poly_interpolate(qM, H, 12);
	poly_interpolate(qC, H, 12);
	poly_interpolate(PI, H, 12);
	
	
	memset(sigma, -1, sizeof(sigma));
	memset(pos, -1, sizeof(pos));
	memset(start, -1, sizeof(start));
	
	for(i=0;i<n;++i) {
		if(pos[va[i]] == -1) start[va[i]] = i;
		else sigma[i] = pos[va[i]];
		pos[va[i]] = i;
	}
	
	for(i=n;i<2*n;++i) {
		if(pos[vb[i-n]] == -1) start[vb[i-n]] = i;
		else sigma[i] = pos[vb[i-n]];
		pos[vb[i-n]] = i;
	}
	
	for(i=2*n;i<3*n;++i) {
		if(pos[vc[i-2*n]] == -1) start[vc[i-2*n]] = i;
		else sigma[i] = pos[vc[i-2*n]];
		pos[vc[i-2*n]] = i;
	}
	
	for(i=0;i<5*round+3;++i) {
		if(start[i] != -1) sigma[start[i]] = pos[i];
	}
	
	fn_set_word(L1[0], 1);
	for(i=1;i<n;++i) {
		fn_set_word(L1[i], 0);
	}
	poly_interpolate(L1, H, 12); // setup
	
	for(i=0;i<n;++i) {
		fn_copy(S_sigma1[i], H[sigma[i]]);
		fn_copy(S_sigma2[i], H[sigma[n+i]]);
		fn_copy(S_sigma3[i], H[sigma[2*n+i]]);
	}
	poly_interpolate(S_sigma1, H, 12); // setup
	poly_interpolate(S_sigma2, H, 12);
	poly_interpolate(S_sigma3, H, 12);
	
	fn_t xl;
	fn_t xr;
	
	fn_from_hex(xr, sage_hex_xr);
	fn_from_hex(xl, sage_hex_xl);
	fn_set_word(wit[0], 0);
	fn_copy(wit[1], xr);
	fn_copy(wit[2], xl);
	
	for(j=0;j<round;++j) {
		fn_copy(wit[5*j+3], constant[j]);
		fn_add(wit[5*j+4], xl, constant[j]);
		fn_sqr(wit[5*j+5], wit[5*j+4]);
		fn_mul(wit[5*j+6], wit[5*j+5], wit[5*j+4]);
		fn_add(wit[5*j+7], wit[5*j+6], xr);
		fn_copy(xr, xl);
		fn_copy(xl, wit[5*j+7]);
	}
	
	for(i=0;i<n;++i) {
		fn_copy(in[i], wit[va[i]]);
		fn_copy(in[n+i], wit[vb[i]]);
		fn_copy(in[2*n+i], wit[vc[i]]);
	}

	plonk_proof(CRS, H, 12, k1, k2, L1, qM, qL, qR, qO, qC, PI,
		sigma, S_sigma1, S_sigma2, S_sigma3, in,
		&A, &B, &C, &Z, &T_lo, &T_mid, &T_hi, &W_delta, &W_delta_w,
		_a, _b, _c, _s_sigma1, _s_sigma2, _t, _r, _z_w);
	
	return 0;
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
	test_plonk_from_sage_simple();
	return 0;
}
