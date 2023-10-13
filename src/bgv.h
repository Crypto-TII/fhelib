#ifndef TIIFHE_BGV_H
#define TIIFHE_BGV_H

#include "config.h"

/***********************
 *   DATA STRUCTURES   *
 ***********************/
typedef struct tiifhe_ciphertext        tiifhe_Ciphertext;
typedef struct tiifhe_ciphertext_switch tiifhe_CiphertextSwitch;
typedef struct tiifhe_delta             tiifhe_Delta;
typedef struct tiifhe_encoding          tiifhe_Encoding;
typedef struct tiifhe_key_public        tiifhe_KeyPublic;
typedef struct tiifhe_key_secret        tiifhe_KeySecret;
typedef struct tiifhe_key_switch        tiifhe_KeySwitch;

struct tiifhe_ciphertext {
	tiifhe_Poly poly[3];
	size_t degree;
	size_t drop;
	int ntt;
};

struct tiifhe_ciphertext_switch {
	tiifhe_Ciphertext ct;
	tiifhe_Poly ext[TIIFHE_OMEGA];
	size_t poly;
};

struct tiifhe_delta {
	tiifhe_Ciphertext ct;
	tiifhe_Digit inv;
	size_t idx;
};

struct tiifhe_encoding {
	tiifhe_Poly value;
	size_t drop;
};

struct tiifhe_key_public {
	tiifhe_Seed seed;
	tiifhe_Poly poly;
};

struct tiifhe_key_secret {
	tiifhe_Poly poly;
};

struct tiifhe_key_switch {
	tiifhe_Seed seed;
	tiifhe_Poly poly[TIIFHE_OMEGA];
};

extern size_t tiifhe_decomp_idx[TIIFHE_OMEGA][TIIFHE_KMAX];
extern size_t tiifhe_decomp_len[TIIFHE_OMEGA];


/**********************
 *   HIGH-LEVEL API   *
 **********************/

/*
 * MEMORY
 *
 * Note that there are no nested allocations for tiifhe_Digit, tiifhe_Poly,
 * and tiifhe_Seed; considering memory layout and efficiency, this assumption
 * should be a reasonable assumption for a ring layer. Memory allocation uses
 * calloc with error checking, thus, allocated memory will be valid and zero.
 */
tiifhe_Ciphertext *tiifhe_bgv_alloc_ct(void);
tiifhe_KeyPublic  *tiifhe_bgv_alloc_pk(void);
tiifhe_KeySecret  *tiifhe_bgv_alloc_sk(void);
tiifhe_KeySwitch  *tiifhe_bgv_alloc_ksw(void);

void tiifhe_bgv_dealloc(void *ptr);
void tiifhe_bgv_dealloc_ct(tiifhe_Ciphertext *ct);
void tiifhe_bgv_dealloc_pk(tiifhe_KeyPublic *pk);
void tiifhe_bgv_dealloc_sk(tiifhe_KeySecret *sk);
void tiifhe_bgv_dealloc_ksw(tiifhe_KeySwitch *ksw);

/* Copy TIIFHE_QPLEN ciphertext structs from ct to rop. */
void tiifhe_bgv_copy(tiifhe_Ciphertext *rop, const tiifhe_Ciphertext *ct);

/*
 * INITIALIZATION
 *
 * These functions initialize the global parameter state and perform
 * the necessary pre-computations, for example, for the NTT. In most cases,
 * the resource cleanup of deinitialization can be deferred to the operating
 * system and no explicit function call is required.
 */
void tiifhe_bgv_deinit(void);
void tiifhe_bgv_init(void);

/*
 * KEY GENERATION
 */
/* Generates a random secret key in the coefficient domain. The default
 * distribution is a centered binomial distribution of width 1, that is
 * a uniform random selection from the set {-1, 0, 0, 1}. */
void tiifhe_bgv_keygen_secret(tiifhe_KeySecret *sk);

/* Generates a public key in the NTT domain for a given secret key. */
void tiifhe_bgv_keygen_public(tiifhe_KeyPublic *pk, const tiifhe_KeySecret *sk);

/* Generates a key switching key for multiplication, that is sk^2. */
void tiifhe_bgv_keygen_switch_mul(tiifhe_KeySwitch *ksw, const tiifhe_KeySecret *sk);

/* Generates a key switching key for rotation, that is rot(sk, steps). */
void tiifhe_bgv_keygen_switch_rot(tiifhe_KeySwitch *ksw, const tiifhe_KeySecret *sk, size_t steps);

/*
 * HOMOMORPHIC OPERATIONS
 */
/* Create a BGV ciphertext given a public key and a message. */
void tiifhe_bgv_encrypt(tiifhe_Ciphertext *ct, tiifhe_KeyPublic *pk, const tiifhe_Message *m);

/* Recover a message given a secret key and ciphertext. */
void tiifhe_bgv_decrypt(tiifhe_Message *m, const tiifhe_KeySecret *sk, tiifhe_Ciphertext *ct);

/* Computes the homomorphic addition rop = op1 + op2. */
void tiifhe_bgv_add(tiifhe_Ciphertext *rop, tiifhe_Ciphertext *op1, tiifhe_Ciphertext *op2);

/* Computes the ciphertext-plaintext addition rop = ct + m. */
void tiifhe_bgv_addc(tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, const tiifhe_Message *m);

/* Computes the homomorphic multiplication rop = op1 * op2. */
void tiifhe_bgv_mul(tiifhe_Ciphertext *rop, tiifhe_Ciphertext *op1, tiifhe_Ciphertext *op2);

/* Computes the homomorphic constant multiplication rop = ct * m. */
void tiifhe_bgv_mulc(tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, const tiifhe_Message *m);

/* Computes the homomorphic rotation rop = rot(ct, steps). */
void tiifhe_bgv_rot(tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, size_t steps);

/*
 * HOMOMORPHIC HOUSEKEEPING
 */
/* Drop the next prime from the ciphertext starting at prime index 0. */
void tiifhe_bgv_modswitch(tiifhe_Ciphertext *ct);

/* Drop a ciphertext down to the level of a comparison ciphertext. */
void tiifhe_bgv_drop(tiifhe_Ciphertext *ct, const tiifhe_Ciphertext *cmp);

/* Perform key switching on a ciphertext after multiplication, that is,
 * re-linearize the ciphertext from three to two active polynomials with ksw
 * as generated via tiifhe_bgv_keygen_switch_mul. */
void tiifhe_bgv_keyswitch_mul(tiifhe_Ciphertext *ct, tiifhe_KeySwitch *ksw);

/* Perform key switching on a ciphertext with two active polynomials, usually
 * after a rotation, with ksw generated via tiifhe_bgv_keygen_switch_rot. */
void tiifhe_bgv_keyswitch_rot(tiifhe_Ciphertext *ct, tiifhe_KeySwitch *ksw);


/*********************
 *   LOW-LEVEL API   *
 *********************/

/* Re-construct a message from a RNS-encoded plaintext. */
void tiifhe_bgv_decode(tiifhe_Message *m, const tiifhe_Encoding *p);

/* Transform the polynomial dsw = 1 or dsw = 2 given a key switching key. */
void tiifhe_bgv_keyswitch(tiifhe_Ciphertext *ct, tiifhe_KeySwitch *ksw, size_t dsw);

/* Perform a base extension for the ciphertext polynomial at index poly. */
void tiifhe_bgv_keyswitch_ext(tiifhe_CiphertextSwitch *csw, tiifhe_Ciphertext *ct, size_t poly);

/* Compute rop = op1 + op2 for the RNS prime at index idx. */
void tiifhe_bgv_idx_add(size_t idx, tiifhe_Ciphertext *rop, tiifhe_Ciphertext *op1, tiifhe_Ciphertext *op2);

/* Compute rop = ct + p for the RNS prime at index idx. */
void tiifhe_bgv_idx_addc(size_t idx, tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, const tiifhe_Encoding *p);

/* Decrypt for the RNS prime at index idx. */
void tiifhe_bgv_idx_decrypt(size_t idx, tiifhe_Encoding *p, const tiifhe_KeySecret *sk, tiifhe_Ciphertext *ct);

/* Encode a message to the RNS prime at index idx. */
void tiifhe_bgv_idx_encode(size_t idx, tiifhe_Encoding *p, const tiifhe_Message *m);

/* Encrypt an encoded message for the RNS prime at index idx. */
void tiifhe_bgv_idx_encrypt(size_t idx, tiifhe_Ciphertext *ct, tiifhe_KeyPublic *pk, const tiifhe_Encoding *p, tiifhe_Seed *seed);

/* Generate a public key for the RNS prime at index idx. Resets the init flag
 * of the public key seed and the error seed after use. Across primes, the seed
 * needs to be consistent for a given encryption. */
void tiifhe_bgv_idx_keygen_public(size_t idx, tiifhe_KeyPublic *pk, const tiifhe_KeySecret *sk, tiifhe_Seed *seed);

/* Generate a key switching key for the RNS prime at index idx. Requires
 * TIIFHE_OMEGA many seeds consistent across all RNS primes; the init flag
 * of each seed is reset after use. Uses sk to hide sk2. */
void tiifhe_bgv_idx_keygen_switch(size_t idx, tiifhe_KeySwitch *ksw, const tiifhe_KeySecret *sk, const tiifhe_KeySecret *sk2, tiifhe_Seed *seed);

/* Calls tiifhe_bgv_idx_keygen_switch with sk2 = sk^2 for the RNS prime
 * at index idx. */
void tiifhe_bgv_idx_keygen_switch_mul(size_t idx, tiifhe_KeySwitch *ksw, const tiifhe_KeySecret *sk, tiifhe_Seed *seed);

/* Calls tiifhe_bgv_idx_keygen_switch with sk2 = rot(sk, steps) for the
 * RNS prime at index idx. */
void tiifhe_bgv_idx_keygen_switch_rot(size_t idx, tiifhe_KeySwitch *ksw, const tiifhe_KeySecret *sk, size_t steps, tiifhe_Seed *seed);

/* Synonym for tiifhe_bgv_idx_modswitch_delta. */
void tiifhe_bgv_idx_keyswitch_delta(size_t idx, tiifhe_Delta *delta, tiifhe_Ciphertext *ct);

/* Perform the dot product for an extended ciphertext with a key switching key
 * for the RNS prime at index idx. */
void tiifhe_bgv_idx_keyswitch_dot(size_t idx, tiifhe_Ciphertext *ct, tiifhe_CiphertextSwitch *csw, tiifhe_KeySwitch *ksw);

/* Base extends delta to TIIFHE_P and calls tiifhe_bgv_idx_modswitch for
 * the RNS prime at index idx. */
void tiifhe_bgv_idx_keyswitch_switch(size_t idx, tiifhe_Ciphertext *ct, tiifhe_Delta *delta);

/* Perform modulus switching given delta for the RNS prime at index idx. */
void tiifhe_bgv_idx_modswitch(size_t idx, tiifhe_Ciphertext *ct, tiifhe_Delta *delta);

/* Compute delta for a given ciphertext for the RNS prime at index idx.
 * Sets drop = 1 for the ciphertext. */
void tiifhe_bgv_idx_modswitch_delta(size_t idx, tiifhe_Delta *delta, tiifhe_Ciphertext *ct);

/* Base extend delta of length len to the RNS prime at index idx and perform
 * modulus switching on the given ciphertext. */
void tiifhe_bgv_idx_modswitch_ext(size_t idx, tiifhe_Ciphertext *ct, tiifhe_Delta *delta, size_t len);

/* Compute rop = op1 * op2 for the RNS prime at index idx. */
void tiifhe_bgv_idx_mul(size_t idx, tiifhe_Ciphertext *rop, tiifhe_Ciphertext *op1, tiifhe_Ciphertext *op2);

/* Compute rop = ct * p for the RNS prime at index idx. */
void tiifhe_bgv_idx_mulc(size_t idx, tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, const tiifhe_Encoding *p);

/* Homomorphically rotate ct by steps and save in rop for the RNS prime at
 * index idx. Requires rop != ct. */
void tiifhe_bgv_idx_rot(size_t idx, tiifhe_Ciphertext *rop, tiifhe_Ciphertext *ct, size_t steps);

/* Homomorphically rotate csw by steps and save in rop for the RNS prime at
 * index idx. Requires rop != csw. */
void tiifhe_bgv_idx_rot_csw(size_t idx, tiifhe_CiphertextSwitch *rop, tiifhe_CiphertextSwitch *csw, size_t steps);

/* Homomorphically rotate csw by steps for the RNS prime at index idx. */
void tiifhe_bgv_idx_rot_csw_inplace(size_t idx, tiifhe_CiphertextSwitch *csw, size_t steps);

/* Homomorphically rotate ct by steps for the RNS prime at index idx. */
void tiifhe_bgv_idx_rot_inplace(size_t idx, tiifhe_Ciphertext *ct, size_t steps);


#endif /* TIIFHE_BGV_H */
