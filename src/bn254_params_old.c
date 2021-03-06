/*
 *  Copyright 2021 The ZKrypt Project. All Rights Reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License"); you may
 *  not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 */

#include "bn254.h"

#if 1
// y^2 = x^3 + 3

const uint32_t BN254_P[] = {
	0xd87cfd47, 0x3c208c16, 0x6871ca8d, 0x97816a91,
	0x8181585d, 0xb85045b6, 0xe131a029, 0x30644e72,
};
const uint32_t BN254_N[] = {
	0xf0000001, 0x43e1f593, 0x79b97091, 0x2833e848,
	0x8181585d, 0xb85045b6, 0xe131a029, 0x30644e72,
};
const uint32_t BN254_P_MU[] = {
	0x9bf90e51, 0xf3aed8a1, 0x7cd4c086, 0xe965e176,
	0x8073013a, 0xb074a586, 0x23a04a7a, 0x4a474626,
};
const uint32_t BN254_N_MU[] = {
	0xe1de9259, 0x20703a6b, 0x9e880ae6, 0x14485200,
	0x80730147, 0xb074a586, 0x23a04a7a, 0x4a474626,
};
const uint32_t BN254_P_INV_NEG[] = {
	0xe4866389, 0x87d20782, 0x1eca6ac9, 0x9ede7d65,
	0x1833da80, 0xd8afcbd0, 0x91888c6b, 0xf57a22b7,
};
const uint32_t BN254_P_MONT_ONE[] = {
	0xc58f0d9d, 0xd35d438d, 0xf5c70b3d, 0xa78eb28,
	0x7879462c, 0x666ea36f, 0x9a07df2f, 0xe0a77c1,
};
const uint32_t BN254_P_MONT_ONE_SQR[] = {
	0x538afa89, 0xf32cfc5b, 0xd44501fb, 0xb5e71911,
	0xa417ff6, 0x47ab1eff, 0xcab8351f, 0x6d89f71,
};
const uint32_t BN254_N_INV_NEG[] = {
	0xe4866389, 0x87d20782, 0x1eca6ac9, 0x9ede7d65,
	0x1833da80, 0xd8afcbd0, 0x91888c6b, 0xf57a22b7,
};
const uint32_t BN254_N_MONT_ONE[] = {
	0x4ffffffb, 0xac96341c, 0x9f60cd29, 0x36fc7695,
	0x7879462e, 0x666ea36f, 0x9a07df2f, 0xe0a77c1,
};
const uint32_t BN254_N_MONT_ONE_SQR[] = {
	0xae216da7, 0x1bb8e645, 0xe35c59e3, 0x53fe3ab1,
	0x53bb8085, 0x8c49833d, 0x7f4e44a5, 0x216d0b1,
};

#ifdef BN254_FP_USE_MONTGOMERY
const point_t BN254_G1 = {
	{ 0xc58f0d9d, 0xd35d438d, 0xf5c70b3d, 0x0a78eb28,
	  0x7879462c, 0x666ea36f, 0x9a07df2f, 0x0e0a77c1 },
	{ 0x8b1e1b3a, 0xa6ba871b, 0xeb8e167b, 0x14f1d651,
	  0xf0f28c58, 0xccdd46de, 0x340fbe5e, 0x1c14ef83 },
	{ 0xc58f0d9d, 0xd35d438d, 0xf5c70b3d, 0x0a78eb28,
	  0x7879462c, 0x666ea36f, 0x9a07df2f, 0x0e0a77c1 },
	0,
};
#else
const point_t BN254_G1 = {
	{ 1,0,0,0,0,0,0,0 },
	{ 2,0,0,0,0,0,0,0 },
	{ 1,0,0,0,0,0,0,0 },
	0,
};
#endif


#else
const uint32_t BN254_P[] = {
	0x00000013, 0xa7000000, 0x00000013, 0x61210000,
	0x00000008, 0xba344d80, 0x40000001, 0x25236482,
};

const uint32_t BN254_N[8] = {
	0x0000000d, 0xa1000000, 0x00000010, 0xff9f8000,
	0x00000007, 0xba344d80, 0x40000001, 0x25236482,
};

// barrett mu(p) = 2^512 // p
const uint32_t BN254_P_MU[] = {
	0x4c735a91, 0x5a5f8d22, 0x3b56f610, 0x24046450,
	0x7bb36c39, 0x7fcedba3, 0x93af3394, 0xe4a64840,
	0x6,
};

// barrett mu(n) = 2^512 // n
const uint32_t BN254_N_MU[9] = {
	0x284df31b, 0x8cbaa5a3, 0x3ad196e5, 0x3d13303c,
	0x7bb36c4b, 0x7fcedba3, 0x93af3394, 0xe4a64840,
	0x06,
};

// montgomery -p^-1 (mod R), R = 2^256
const uint32_t BN254_P_INV_NEG[] = {
	0xd79435e5, 0x08435e50, 0x1104f6c8, 0x6e371ba8,
	0xc45b843c, 0x92022379, 0xba60808c, 0xb65373cc,
};

// montgomery mont(1) (mod p) = 1*R (mod p) = R (mod p), R = 2^256
const uint32_t BN254_P_MONT_ONE[] = {
	0xffffff8e, 0x15ffffff, 0xffffff8a, 0xb939ffff,
	0xffffffcd, 0xa2c62eff, 0x7ffffff5, 0x212ba4f2,
};

// montgomery R^2 (mod p), R = 2^256
const uint32_t BN254_P_MONT_ONE_SQR[] = {
	0x5370473d, 0xb3e88674, 0x8c1cc3f1, 0x55efbf6e,
	0x7f86954f, 0x281e3a1b, 0xf6403a3d, 0x1b0a32fd,
};

// montgomery -n^-1 (mod R), R = 2^256
const uint32_t BN254_N_INV_NEG[] = {
	0x3b13b13b, 0xea3b13b1, 0x731fcf86, 0x4b35ee94,
	0x1df24016, 0x42217fc3, 0xa9bba566, 0xf3d5266b,
};

// montgomery mont(1) (mod n) = 1*R (mod n) = R (mod n), R = 2^256
const uint32_t BN254_N_MONT_ONE[] = {
	0xffffffb2, 0x39ffffff, 0xffffff9c, 0x0242ffff,
	0xffffffd0, 0xa2c62eff, 0x7ffffff5, 0x212ba4f2,
};

// montgomery R^2 (mod n), R = 2^256
const uint32_t BN254_N_MONT_ONE_SQR[] = {
	0xf40aa7a1, 0xdf8596b6, 0xe2231ec3, 0xe0885092,
	0x575d5a78, 0xc300765b, 0x325f9035, 0x24e8b3bc,
};

#ifdef BN254_FP_USE_MONTGOMERY
const point_t BN254_G1 = {
	{ 0x00000085, 0x91000000, 0x00000089, 0xa7e70000,
	  0x0000003a, 0x176e1e80, 0xc000000c, 0x03f7bf8f },
	{ 0xffffff8e, 0x15ffffff, 0xffffff8a, 0xb939ffff,
	  0xffffffcd, 0xa2c62eff, 0x7ffffff5, 0x212ba4f2 },
	{ 0xffffff8e, 0x15ffffff, 0xffffff8a, 0xb939ffff,
	  0xffffffcd, 0xa2c62eff, 0x7ffffff5, 0x212ba4f2 },
	0,
};
#else
const point_t BN254_G1 = {
	{ 0x00000012, 0xa7000000, 0x00000013, 0x61210000,
	  0x00000008, 0xba344d80, 0x40000001, 0x25236482 },
	{ 1,0,0,0,0,0,0,0 },
	{ 1,0,0,0,0,0,0,0 },
	0,
};
#endi

#endif
